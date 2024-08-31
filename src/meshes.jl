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
        if _f isa AbstractVector{<: MultiFace}
            if facetype isa MultiFace
                # drop faces that facetype doesn't include
                _f = simplify_faces(facetype, _f)
            else
                # drop faces for vertex attributes that aren't given
                names = (:position, keys(vertex_attributes)...)
                _f2 = simplify_faces(names, _f)
                
                # and remap to a simple face type so that decompose can handle the rest
                _f, mappings = merge_vertex_indices(_f2)
                positions = positions[mappings[1]]
                vertex_attributes = NamedTuple(
                    (Pair(names[i], vertex_attributes[i-1][mappings[i]]) for i in 2:length(mappings))    
                )
            end
        end
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

function triangle_mesh(primitive::Union{AbstractGeometry{N}, AbstractVector{<: Point{N}}})::SimpleTriangleMesh{N} where {N}
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
        primitive, uv = decompose_uv(primitive), normal = decompose_normals(primitive),
        pointtype=Point{N,Float32}, facetype=GLTriangleFace)
end

function normal_mesh(points::AbstractVector{<:Point},
                     faces::AbstractVector{<:AbstractFace})
    _points = decompose(Point3f, points)
    _faces = decompose(GLTriangleFace, faces)
    return Mesh(_faces, position = _points, normal = normals(_points, _faces))
end

function normal_mesh(primitive::AbstractGeometry{N}) where {N}
    return mesh(
        primitive, normal = decompose_normals(primitive), 
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

# TODO: Is this ok as "public" function?
# MultiFace(f1, f2, f3) + (o1, o2, o3) = MultiFace(f1 + o1, f2 + o2, f3 + o3)
function Base.:+(f::MultiFace{N, T, FT, Names, M}, o::NTuple{M, T}) where {N, T, FT, Names, M}
    return MultiFace{Names}(ntuple(m -> f[m] + o[m], M))
end

function Base.merge(meshes::AbstractVector{<:Mesh})
    return if isempty(meshes)
        return Mesh(Point3f[], GLTriangleFace[])
    elseif length(meshes) == 1
        return meshes[1]
    else

        m1 = meshes[1]

        # Check that all meshes use the same VertexAttributes
        # Could also do this via typing the function, but maybe error is nice?
        names = propertynames(m1.vertex_attributes)
        idx = findfirst(m -> propertynames(m.vertex_attributes) != names, meshes)
        if idx !== nothing
            error(
                "Cannot merge meshes with different vertex attributes. " * 
                "First missmatch between meshes[1] with $names and " * 
                "meshes[$idx] with $(propertynames(meshes[idx]))."
            )
        end

        # We can't merge MultiFace with standard faces because MutliFace allows
        # desynchronizes vertex indices that normal faces assume synchronized.
        is_multi = facetype(m1) <: MultiFace

        if all(m -> is_multi == (facetype(m) <: MultiFace), meshes)

            # All the same kind of face, can just merge

            new_attribs = NamedTuple{names}(map(names) do name
                return mapreduce(m -> getproperty(m, name), vcat, meshes)
            end)
            fs = reduce(vcat, faces.(meshes))

            # TODO: is the type difference in offset bad?
            idx = length(faces(m1))
            offset = is_multi ? length.(values(vertex_attributes(m1))) : length(coordinates(m1))
            views = isempty(m1.views) ? [1:idx] : copy(m1.views)
            
            for mesh in Iterators.drop(meshes, 1)
                # update face indices
                N = length(faces(mesh))
                for i = idx .+ (1:N)
                    fs[i] = fs[i] + offset
                end

                # add views
                if isempty(mesh.views)
                    push!(views, idx+1 : idx+N)
                else
                    append!(views, (view + idx for view in mesh.views))
                end
    
                idx += N
                if is_multi
                    offset = offset .+ length.(values(vertex_attributes(mesh)))
                else
                    offset += length(coordinates(mesh))
                end
            end

            return Mesh(new_attribs, fs, views)

        else

            # TODO: We can probably simplify this to `merge(merge_vertex_indices.(meshes))`
            #       but need to check performance


            # Varying ace types, need to convert MultiFace
            new_attribs = NamedTuple{names}(similar.(values(vertex_attributes(m1)), 0))
            FT = facetype(m1) <: MultiFace ? eltype(facetype(m1)) : facetype(m1)
            remapped_faces = []
            new_views = UnitRange{Int}[]
            vertex_index_counter = eltype(FT)(1)

            for mesh in meshes
                # convert MultiFace mesh to normal faces, synchronizing vertex indices
                attribs, fs, views = merge_vertex_indices(
                    vertex_attributes(mesh), faces(mesh), mesh.views, vertex_index_counter)
                
                # increment first vertex index used by faces of the next iteration
                vertex_index_counter += length(attribs[1])
                
                # update merged data
                for name in names
                    append!(new_attribs[name], attribs[name])
                end

                push!(remapped_faces, fs)

                if isempty(views)
                    push!(new_views, 1:length(fs))
                else 
                    append!(new_views, views)
                end
            end

            # We did MultiFace -> normal face, now equalize normal face types
            new_faces = reduce(vcat, remapped_faces)

            return Mesh(new_attribs, new_faces, new_views)
        end

    end
end

# TODO: Probably not our problem
# function Base.merge(meshes::AbstractVector{T}) where T <: MetaMesh
#     isempty(meshes) && return T(Point3f[], GLTriangleFace[])
#     big_mesh = merge(map(Mesh, meshes))
#     big_meta = deepcopy(meta(meshes[1]))
#     for mesh in Iterators.drop(meshes, 1)
#         mm = meta(mesh)
#         for (k, v) in pairs(mm)
#             append!(big_meta[k], v)
#         end
#     end
#     return MetaMesh(big_mesh, big_meta)
# end

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
        faces::AbstractVector{<: FT},
        views::Vector{UnitRange{Int}},
        vertex_index_counter = nothing
    ) where {Names, FT <: AbstractFace}

    if FT <: MultiFace
        error(
            "Failed to call correct method. This likely happened because vertex " *
            "attributes names $Names do not match face name $(propertynames(first(faces)))."
        )
    end

    return attribs, faces, views
end

function merge_vertex_indices(
        attribs::NamedTuple{Names}, 
        faces::AbstractVector{<: MultiFace{N, T, FT, Names}},
        views::Vector{UnitRange{Int}},
        vertex_index_counter = T(1) # TODO: test 0 vs 1 base
    ) where {Names, N, T, FT}

    # Note: typing checks for matching Names

    if isempty(views)
        new_faces, vertex_map = merge_vertex_indices(faces)
        new_attribs = ntuple(n -> attribs[n][vertex_map[n]], length(Names))
        return NamedTuple{Names}(new_attribs), new_faces, views
    end

    new_attribs = NamedTuple((Pair(k, similar(v, 0)) for (k, v) in pairs(attribs)))
    new_faces = similar(faces, FT, 0)
    new_views = UnitRange{Int}[]

    # TODO: this depends on T in Face (1 based -> 1, 0 based -> 0)
    for idxs in views
        # Generate new face from current view, with the first vertex_index 
        # corresponding to the first vertex attribute added in this iteration
        face_view = view(faces, idxs)
        new_faces_in_view, vertex_map = merge_vertex_indices(face_view, vertex_index_counter)
        vertex_index_counter += length(vertex_map)

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
        faces::AbstractVector{FT},
        vertex_index_counter = T(1)
    ) where {N, T, FT <: AbstractFace{N, T}}

    @assert !(FT <: MultiFace) "Dispatch failed?"
    
    N_vert = mapreduce(f -> max(f), max, faces)
    return faces, (1:N_vert) .+ vertex_index_counter
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
            vertex = ntuple(n -> multi_face[n][i], N_Attrib)

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
