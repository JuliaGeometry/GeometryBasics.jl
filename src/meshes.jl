
"""
Mesh <: AbstractMesh{Element}
The concrete AbstractMesh type.
"""
struct Mesh{Dim, T<:Number, V<:AbstractVector{Point{Dim, T}}, C <: AbstractVector{<: AbstractFace}} <: AbstractMesh{Dim, T}
    vertices::V
    connectivity::C
end
coordinates(mesh::Mesh) = mesh.vertices
faces(mesh::Mesh) = mesh.connectivity
Base.getindex(mesh::Mesh, i::Integer) = mesh.vertices[mesh.connectivity[i]]
Base.length(mesh::Mesh) = length(mesh.connectivity)
Base.:(==)(a::Mesh, b::Mesh) = coordinates(a) == coordinates(b) && faces(a) == faces(b)

function Base.iterate(mesh::Mesh, i=1)
    return i - 1 < length(mesh) ? (mesh[i], i + 1) : nothing
end

function Mesh(points::AbstractVector{Point{Dim, T}},
          faces::AbstractVector{<:AbstractFace}) where {Dim, T}
    return Mesh{Dim, T, }(points, faces)
end

function Mesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:Integer},
          facetype=TriangleFace, skip=1)
    return Mesh(points, connect(faces, facetype, skip))
end

struct MetaMesh{Dim, T, M <: AbstractMesh{Dim, T}, Names, Types} <: AbstractMesh{Dim, T}
    mesh::M
    meta::NamedTuple{Names, Types}
    function MetaMesh(mesh::AbstractMesh{Dim, T}, meta::NamedTuple{Names, Types}) where {Dim, T, Names, Types}
        new{Dim, T, typeof(mesh), Names, Types}(mesh, meta)
    end
end

function MetaMesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractFace}; meta...)
    MetaMesh(Mesh(points, faces), values(meta))
end

function MetaMesh(m::AbstractMesh; meta...)
    MetaMesh(m, values(meta))
end

@inline Base.hasproperty(mesh::MetaMesh, field::Symbol) = hasproperty(getfield(mesh, :meta), field)
@inline Base.getproperty(mesh::MetaMesh, field::Symbol) = getproperty(getfield(mesh, :meta), field)
@inline Base.propertynames(mesh::MetaMesh) = propertynames(getfield(mesh, :meta))

coordinates(mesh::MetaMesh) = coordinates(Mesh(mesh))
faces(mesh::MetaMesh) = faces(Mesh(mesh))
normals(mesh::MetaMesh) = hasproperty(mesh, :normals) ? mesh.normals : nothing
texturecoordinates(mesh::MetaMesh) = hasproperty(mesh, :uv) ? mesh.uv : nothing

meta(mesh::MetaMesh) = getfield(mesh, :meta)
Mesh(mesh::MetaMesh) = getfield(mesh, :mesh)

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

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)

Create a mesh from a polygon given as a vector of points, using triangulation.
"""
function mesh(polygon::AbstractVector{P}; pointtype=P) where {P<:Point{2}}
    return mesh(Polygon(polygon); pointtype=pointtype)
end

function triangle_mesh(primitive::AbstractGeometry{N}) where {N}
    return mesh(primitive; pointtype=Point{N, Float32})
end

function triangle_mesh(primitive::AbstractVector{<: Point})
    return mesh(Polygon(primitive))
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
        m1 = meshes[1]
        ps = copy(coordinates(m1))
        fs = copy(faces(m1))
        for mesh in Iterators.drop(meshes, 1)
            append!(fs, map(f -> f .+ length(ps), faces(mesh)))
            append!(ps, coordinates(mesh))
        end
        return Mesh(ps, fs)
    end
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
