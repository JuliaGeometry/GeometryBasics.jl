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
function mesh(primitive::AbstractGeometry; pointtype=Point, facetype=GLTriangleFace, vertex_attributes...)
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
    return Mesh(positions, f; vertex_attributes...)
end

const SimpleMesh{N, T, FT} = Mesh{N, T, FT, (:position,), Tuple{Vector{Point{N, T}}}, Vector{FT}}
const SimpleTriangleMesh{N} = SimpleMesh{N, Float32, GLTriangleFace}

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

function triangle_mesh(primitive::Union{AbstractGeometry{N}, AbstractVector{<: Point{N}}}; nvertices = nothing)::SimpleTriangleMesh{N} where {N}
    if nvertices !== nothing
        @warn("nvertices argument deprecated. Wrap primitive in `Tesselation(primitive, nvertices)`")
        primitive = Tesselation(primitive, nvertices)
    end
    return mesh(primitive; pointtype=Point{N,Float32}, facetype=GLTriangleFace)
end

function uv_mesh(primitive::AbstractGeometry{N,T}) where {N,T}
    return mesh(primitive, uv = decompose_uv(primitive), pointtype=Point{N,Float32}, facetype=GLTriangleFace)
end

function uv_normal_mesh(primitive::AbstractGeometry{N}) where {N}
    return mesh(
        primitive, uv = decompose_uv(primitive), normals = decompose_normals(primitive),
        pointtype=Point{N,Float32}, facetype=GLTriangleFace)
end

function normal_mesh(points::AbstractVector{<:Point},
                     faces::AbstractVector{<:AbstractFace})
    _points = decompose(Point3f, points)
    _faces = decompose(GLTriangleFace, faces)
    return Mesh(_faces, position = _points, normals = normals(_points, _faces))
end

function normal_mesh(primitive::AbstractGeometry{N}) where {N}
    return mesh(
        primitive, normals = decompose_normals(primitive), 
        pointtype=Point{N,Float32}, facetype=GLTriangleFace)
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

# TODO: naming
# synchronize_vertex_attributes
# merge_vertex_(attribute)_indices
# convert(Face, MultiFace)
# ...
function merge_vertex_indices(mesh)
    attribs, fs, views = merge_vertex_indices(
        vertex_attributes(mesh), faces(mesh), mesh.views)

    return Mesh(attribs, fs, views)
end

function merge_vertex_indices(
        attribs::NamedTuple{Names}, 
        faces::AbstractVector{<: MultiFace{N, T, FT, Names}},
        views::Vector{UnitRange}
    ) where {Names, N, T, FT}

    # Note: typing checks for matching Names

    if isempty(views)
        new_faces, vertex_map = merge_vertex_indices(faces)
        new_attribs = ntuple(n -> attribs[n][vertex_map[n]], length(Names))
        return NamedTuple{Names}(new_attribs), new_faces, views
    end

    new_attribs = NamedTuple((Pair(k, similar(v, 0)) for (k, v) in pairs(attribs)))
    new_faces = similar(faces, FT, 0)
    new_views = UnitRange[]

    for idxs in views
        # TODO: this depends on T in Face (1 based -> 1, 0 based -> 0)
        vertex_index_counter = T(length(new_attribs[1]) + 1)

        # Generate new face from current view, with the first vertex_index 
        # corresponding to the first vertex attribute added in this iteration
        face_view = view(faces, idxs)
        new_faces_in_view, vertex_map = merge_vertex_indices(face_view, vertex_index_counter)
        
        # update vertex attributes
        for (name, indices) in pairs(vertex_map)
            append!(new_attribs[name], view(attribs[name], indices))
        end

        # add new faces and new view
        start = length(new_faces) + 1
        append!(new_faces, new_faces_in_view)
        append!(new_views, start:length(new_faces))
    end

    return new_attribs, new_faces, new_views
end

function merge_vertex_indices(
        faces::AbstractVector{<: MultiFace{N, T, FT, Names, N_Attrib}},
        vertex_index_counter = T(1)
    ) where {N, T, FT <: AbstractFace{N, T}, Names, N_Attrib}

    N_faces = length(faces)

    # maps a combination of old indices in MultiFace to a new vertex_index
    vertex_index_map = Dict{NTuple{N_Attrib, T}, T}()

    # Faces after conversion
    new_faces = sizehint!(FT[], N_faces)

    # indices that remap attributes
    attribute_indices = ntuple(n -> sizehint!(UInt32[], N_faces), N_Attrib)

    # keep track of the remmaped indices for one vertex so we don't have to 
    # query the dict twice
    temp = zeros(N)

    for multi_face in faces

        for i in 1:N
            # get the i-th set of vertex indices from multi_face, i.e.
            # (multi_face.position_index[i], multi_face.normal_index[i], ...)
            vertex = ntuple(n -> multi_face.faces[n][i], N_Attrib)

            # if the vertex exists, get it's index
            # otherwise register it with the next available vertex index
            temp[i] = get!(vertex_index_map, vertex) do 
                vertex_index_counter += 1
                push!.(attribute_indices, vertex)
                return vertex_index_counter - 1
            end
        end

        # generate new face
        push!(new_faces, FT(temp))
    end

    # in case we are reserving more than needed
    sizehint!(new_faces, length(new_faces))

    return new_faces, attribute_indices
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

# TODO:
add_meta(m, kw...) = error("TODO")
pop_meta(m, kw...) = error("TODO")
# function add_meta(mesh::MetaMesh; kw...)
#     return MetaMesh(Mesh(mesh), (; meta(mesh)..., kw...))
# end

# function add_meta(mesh::Mesh; kw...)
#     return MetaMesh(mesh, (; meta(mesh)..., kw...))
# end

# # I didn't find a simple way to remove a field from a namedtuple in a type stable way without
# # a generated function..
# @generated function pop(nt::NamedTuple{Names, Values}, ::Val{name}) where {Names, Values, name}
#     if !(name in Names)
#         return :(throw(Base.KeyError($(QuoteNode(name)))))
#     else
#         names = filter(x-> x !== name, Names)
#         nt = map(names) do name
#             :($name = nt.$(name))
#         end
#         return :((; $(nt...)), nt.$(name))
#     end
# end

# function pop_meta(mesh::MetaMesh, name::Symbol)
#     new_meta, value = pop(meta(mesh), Val(name))
#     return MetaMesh(mesh, new_meta), value
# end

function Base.show(io::IO, ::MIME"text/plain", mesh::Mesh{N, T, FT}) where {N, T, FT}
    println(io, "Mesh{$N, $T, $FT}")
    println(io, "    faces: ", length(faces(mesh)))
    for (name, attrib) in pairs(vertex_attributes(mesh))
        println(io, "    vertex $(name): ", length(attrib))
    end
end

function Base.show(io::IO, ::Mesh{N, T, FT}) where {N, T, FT}
    print(io, "Mesh{$N, $T, $(FT)}(...)")
end

function Base.show(io::IO, ::MIME"text/plain", mesh::MetaMesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    println(io, "MetaMesh{$N, $T, $(FT)}")
    println(io, "    faces: ", length(faces(mesh)))
    for (name, attrib) in pairs(vertex_attributes(mesh))
        println(io, "    vertex $(name): ", length(attrib))
    end
    println(io, "    meta: ", keys(mesh.meta))
end

function Base.show(io::IO, mesh::MetaMesh{N, T}) where {N, T}
    FT = eltype(faces(mesh))
    println(io, "MetaMesh{$N, $T, $(FT)}($(join(keys(meta(mesh)), ", ")))")
end
