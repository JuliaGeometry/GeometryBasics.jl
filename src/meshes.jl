"""
    mesh(primitive::GeometryPrimitive;
         pointtype=Point, facetype=GLTriangle,
         uv=nothing, normaltype=nothing)

Creates a mesh from `primitive`.

Uses the element types from the keyword arguments to create the attributes.
The attributes that have their type set to nothing are not added to the mesh.
Note, that this can be an `Int` or `Tuple{Int, Int}``, when the primitive is grid based.
It also only losely correlates to the number of vertices, depending on the algorithm used.
#TODO: find a better number here!
"""
function mesh(primitive::AbstractGeometry; pointtype=Point, facetype=GLTriangleFace)
    positions = decompose(pointtype, primitive)
    _f = faces(primitive)
    # If faces returns nothing for primitive, we try to triangulate!
    if isnothing(_f)
        if eltype(positions) <: Point2
            # triangulation.jl
            f = decompose(facetype, positions)
        else
            error("No triangulation for $(typeof(primitive))")
        end
    else
        f = decompose(facetype, _f)
    end
    return Mesh(positions, f)
end

const SimpleMesh{N, T, FT} = Mesh{N, T, Vector{Point{N, T}}, Vector{FT}}
const TriangleMesh{N} = SimpleMesh{N, Float32, GLTriangleFace}

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)

Create a mesh from a polygon given as a vector of points, using triangulation.
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace) where {P<:Point{2}}
    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype)
end

function triangle_mesh(primitive::Union{AbstractGeometry{N}, AbstractVector{<: Point{N}}})::TriangleMesh{N} where {N}
    return mesh(primitive; pointtype=Point{N, Float32})
end


pointtype(x::Mesh) = eltype(decompose(Point, x))
facetype(x::Mesh) = eltype(faces(x))

function triangle_mesh(primitive::Mesh{N}) where {N}
    # already target type:
    if pointtype(primitive) === Point{N,Float32} && GLTriangleFace === facetype(primitive)
        return primitive
    else
        return mesh(primitive; pointtype=Point{N,Float32}, facetype=GLTriangleFace)
    end
end

function uv_mesh(primitive::AbstractGeometry{N,T}) where {N,T}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (uv=decompose_uv(primitive),))
end

function uv_normal_mesh(primitive::AbstractGeometry{N}) where {N}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (uv=decompose_uv(primitive), normals=decompose_normals(m)))
end

function normal_mesh(points::AbstractVector{<:Point},
                     faces::AbstractVector{<:AbstractFace})
    _points = decompose(Point3f, points)
    _faces = decompose(GLTriangleFace, faces)
    return MetaMesh(Mesh(_points, _faces), (normals=normals(_points, _faces),))
end

function normal_mesh(primitive::AbstractGeometry{N}) where {N}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (normals=decompose_normals(m),))
end

"""
    volume(triangle)

Calculate the signed volume of one tetrahedron. Be sure the orientation of your
surface is right.
"""
function volume(triangle::Triangle)
    v1, v2, v3 = triangle
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ (v2 × v3)) / 6
end

"""
    volume(mesh)

Calculate the signed volume of all tetrahedra. Be sure the orientation of your
surface is right.
"""
function volume(mesh::Mesh)
    return sum(volume, mesh)
end

function Base.merge(meshes::AbstractVector{<:Mesh})
    return if isempty(meshes)
        return Mesh(Point3f[], GLTriangleFace[])
    elseif length(meshes) == 1
        return meshes[1]
    else
        ps = reduce(vcat, coordinates.(meshes))
        fs = reduce(vcat, faces.(meshes))
        idx = length(faces(meshes[1]))
        offset = length(coordinates(meshes[1]))
        for mesh in Iterators.drop(meshes, 1)
            N = length(faces(mesh))
            for i = idx .+ (1:N)
                fs[i] = fs[i] .+ offset
            end
            idx += N
            offset += length(coordinates(mesh))
        end
        return Mesh(ps, fs)
    end
end

function Base.merge(meshes::AbstractVector{T}) where T <: MetaMesh
    isempty(meshes) && return T(Point3f[], GLTriangleFace[])
    big_mesh = merge(map(Mesh, meshes))
    big_meta = deepcopy(meta(meshes[1]))
    for mesh in Iterators.drop(meshes, 1)
        mm = meta(mesh)
        for (k, v) in pairs(mm)
            append!(big_meta[k], v)
        end
    end
    return MetaMesh(big_mesh, big_meta)
end


function map_coordinates(f, mesh::Mesh)
    result = copy(mesh)
    map_coordinates!(f, result)
    return result
end

function map_coordinates(f, mesh::MetaMesh)
    result = copy(Mesh(mesh))
    map_coordinates!(f, result)
    return MetaMesh(result, meta(mesh))
end

function map_coordinates!(f, mesh::AbstractMesh)
    points = coordinates(mesh)
    map!(f, points, points)
    return mesh
end

function add_meta(mesh::MetaMesh; kw...)
    return MetaMesh(Mesh(mesh), (; meta(mesh)..., kw...))
end

function add_meta(mesh::Mesh; kw...)
    return MetaMesh(mesh, (; meta(mesh)..., kw...))
end

# I didn't find a simple way to remove a field from a namedtuple in a type stable way without
# a generated function..
@generated function pop(nt::NamedTuple{Names, Values}, ::Val{name}) where {Names, Values, name}
    if !(name in Names)
        return :(throw(Base.KeyError($(QuoteNode(name)))))
    else
        names = filter(x-> x !== name, Names)
        nt = map(names) do name
            :($name = nt.$(name))
        end
        return :((; $(nt...)), nt.$(name))
    end
end

function pop_meta(mesh::MetaMesh, name::Symbol)
    new_meta, value = pop(meta(mesh), Val(name))
    return MetaMesh(mesh, new_meta), value
end


function Base.get(f, mesh::MetaMesh, key::Symbol)
    hasproperty(mesh, key) && return getproperty(mesh, key)
    return f()
end

function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    println(io, "Mesh{$N, $T, $(FT)}")
    println(io, "    vertices: ", length(coordinates(mesh)))
    println(io, "    faces: ", length(faces(mesh)), " $(FT)")
end

function Base.show(io::IO, mesh::Mesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    print(io, "Mesh{$N, $T, $(FT)}(...)")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::MetaMesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    println(io, "MetaMesh{$N, $T, $(FT)}")
    println(io, "    vertices: ", length(coordinates(mesh)))
    println(io, "    faces: ", length(faces(mesh)), " $(FT)")
    for (k, v) in pairs(meta(mesh))
        println(io, "    ", k, ": ", length(v), " $(eltype(v))")
    end
end

function Base.show(io::IO, mesh::MetaMesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    println(io, "MetaMesh{$N, $T, $(FT)}($(join(keys(meta(mesh)), ", ")))")
end
