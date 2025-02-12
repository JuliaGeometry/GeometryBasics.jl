"""
    mesh(primitive::GeometryPrimitive[; pointtype = Point, facetype = GLTriangleFace, vertex_attributes...])

Creates a mesh from a given `primitive` with the given `pointtype` and `facetype`.

This method only generates positions and faces from the primitive. Additional
vertex attributes like normals and texture coordinates can be given as extra
keyword arguments.

Note that vertex attributes that are `nothing` get removed before creating a mesh.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function mesh(primitive::AbstractGeometry; pointtype=Point, facetype=GLTriangleFace, vertex_attributes...)
    positions = decompose(pointtype, primitive)

    # TODO: consider not allowing FaceView here?
    if positions isa FaceView
        positions = positions.data
        _fs = positions.faces
        isnothing(faces(primitive)) || @error("A primitive should not define `faces` and use a FaceView for `coordinates()`. Using faces from FaceView.")
    else
        # This tries `faces(prim)` first, then triangulation with the natural
        # position type of the primitive.
        _fs = decompose(facetype, primitive)
    end

    # If faces returns nothing for primitive, we try to triangulate!
    if isnothing(_fs)
        if eltype(positions) <: Point2
            # try triangulation with the converted positions as a last attempt
            fs = decompose(facetype, positions)
        else
            error("No triangulation for $(typeof(primitive))")
        end
    else
        fs = _fs
    end

    return mesh(positions, collect(fs); facetype = facetype, vertex_attributes...)
end

"""
    mesh(positions, faces[, facetype = GLTriangleFace, vertex_attributes...])

Creates a mesh from the given positions and faces. Other vertex attributes like
normals and texture coordinates can be added as keyword arguments.

Note that vertex attributes that are `nothing` get removed before creating a mesh.

See also:[`normal_mesh`](@ref)
"""
function mesh(
        positions::AbstractVector{<:Point},
        faces::AbstractVector{FT};
        facetype=GLTriangleFace, vertex_attributes...
    ) where {FT <: AbstractFace}

    fs = decompose(facetype, faces)

    names = keys(vertex_attributes)
    valid_names = filter(name -> !isnothing(vertex_attributes[name]), names)
    vals = convert_facetype.(Ref(facetype), getindex.(Ref(vertex_attributes), valid_names))
    va = NamedTuple{valid_names}(vals)

    return Mesh(positions, fs; va...)
end

"""
    mesh(mesh::Mesh[; pointtype = Point, facetype = GLTriangleFace, vertex_attributes...]

Recreates the given `mesh` with the given `pointtype`, `facetype` and vertex
attributes. If the new mesh would match the old mesh, the old mesh is returned instead.

Note that vertex attributes that are `nothing` get removed before creating a mesh.
"""
function mesh(
        mesh::Mesh{D, T, FT}; pointtype = Point{D, Float32},
        facetype::Type{<: AbstractFace} = GLTriangleFace,
        attributes...
    ) where {D, T, FT <: AbstractFace}

    names = keys(attributes)
    valid_names = filter(name -> !isnothing(attributes[name]), names)

    if isempty(valid_names) && (GeometryBasics.pointtype(mesh) == pointtype) && (FT == facetype)
        return mesh
    else
        vals = getindex.(Ref(attributes), valid_names)
        va = NamedTuple{valid_names}(vals)

        # add vertex attributes
        va = merge(vertex_attributes(mesh), va)
        # convert position attribute and facetypes in FaceViews
        va = NamedTuple{keys(va)}(map(keys(va)) do name
            val = name == :position ? decompose(pointtype, va[:position]) : va[name]
            return convert_facetype(facetype, val)
        end)

        # update main face type
        f, views = decompose(facetype, faces(mesh), mesh.views)
        return Mesh(va, f, views)
    end
end

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace)

Create a mesh from a polygon given as a vector of points, using triangulation.
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace) where {P<:Point{2}}
    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype)
end

"""
    triangle_mesh(primitive::GeometryPrimitive[; pointtype = Point, facetype = GLTriangleFace])

Creates a simple triangle mesh from a given `primitive` with the given `pointtype`
and `facetype`.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function triangle_mesh(
        primitive::Union{AbstractGeometry{N}, AbstractVector{<: Point{N}}};
        pointtype = Point{N, Float32}, facetype = GLTriangleFace
    ) where {N}
    return mesh(primitive; pointtype = pointtype, facetype = facetype)
end

pointtype(::Mesh{D, T}) where {D, T} = Point{D, T}
facetype(::Mesh{D, T, FT}) where {D, T, FT} = FT

"""
    uv_mesh(primitive::GeometryPrimitive{N}[; pointtype = Point{N, Float32}, facetype = GLTriangleFace, uvtype = Vec2f])

Creates a triangle mesh with texture coordinates from a given `primitive`. The
`pointtype`, `facetype` and `uvtype` are set by the corresponding keyword arguments.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function uv_mesh(
        primitive::AbstractGeometry{N}; pointtype = Point{N, Float32},
        uvtype = Vec2f, facetype = GLTriangleFace
    ) where {N}

    return mesh(
        primitive, uv = decompose_uv(uvtype, primitive),
        pointtype = pointtype, facetype = facetype
    )
end

"""
    uv_normal_mesh(primitive::GeometryPrimitive{N}[; pointtype = Point{N, Float32}, facetype = GLTriangleFace, uvtype = Vec2f, normaltype = Vec3f])

Creates a triangle mesh with texture coordinates and normals from a given
`primitive`. The `pointtype`, `facetype` and `uvtype` and `normaltype` are set
by the corresponding keyword arguments.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function uv_normal_mesh(
        primitive::AbstractGeometry{N}; pointtype = Point{N, Float32},
        uvtype = Vec2f, normaltype = Vec3f, facetype = GLTriangleFace
    ) where {N}

    return mesh(
        primitive, uv = decompose_uv(uvtype, primitive),
        normal = decompose_normals(normaltype, primitive),
        pointtype = pointtype, facetype = facetype
    )
end

"""
    uv_normal_mesh(primitive::GeometryPrimitive{N}[; pointtype = Point{N, Float32}, facetype = GLTriangleFace, uvtype = Vec2f, normaltype = Vec3f])

Creates a triangle mesh with texture coordinates and normals from a given
`primitive`. The `pointtype`, `facetype` and `uvtype` and `normaltype` are set
by the corresponding keyword arguments.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function normal_mesh(
        points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractFace};
        pointtype = Point3f, normaltype = Vec3f, facetype = GLTriangleFace
    )
    _points = decompose(pointtype, points)
    _faces = decompose(facetype, faces)
    return Mesh(_faces, position = _points, normal = normals(_points, _faces, normaltype))
end

"""
    normal_mesh(primitive::GeometryPrimitive{N}[; pointtype = Point{N, Float32}, facetype = GLTriangleFace, normaltype = Vec3f])

Creates a triangle mesh with normals from a given `primitive`. The `pointtype`, `facetype` and `uvtype` and `normaltype` are set
by the corresponding keyword arguments.

See also: [`triangle_mesh`](@ref), [`normal_mesh`](@ref), [`uv_mesh`](@ref), [`uv_normal_mesh`](@ref)
"""
function normal_mesh(
        primitive::AbstractGeometry{N}; pointtype = Point{N, Float32},
        normaltype = Vec3f, facetype = GLTriangleFace
    ) where {N}

    return mesh(
        primitive, normal = decompose_normals(normaltype, primitive),
        pointtype = pointtype, facetype = facetype)
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


"""
    merge(meshes::AbstractVector{Mesh})

Generates a new mesh containing all the data of the individual meshes.

If all meshes are consistent in their use of FaceViews they will be preserved.
Otherwise all of them will be converted with `expand_faceviews(mesh)`.

This function will generate `views` in the new mesh which correspond to the
inputs of this function.
"""
function Base.merge(meshes::AbstractVector{<:Mesh})
    return if isempty(meshes)
        return Mesh(Point3f[], GLTriangleFace[])

    elseif length(meshes) == 1
        return meshes[1]

    else
        m1 = meshes[1]

        # Check that all meshes use the same VertexAttributes
        # Could also do this via typing the function, but maybe error is nice?
        names = keys(vertex_attributes(m1))
        if !all(m -> keys(vertex_attributes(m)) == names, meshes)
            idx = findfirst(m -> keys(vertex_attributes(m)) != names, meshes)
            error(
                "Cannot merge meshes with different vertex attributes. " *
                "First missmatch between meshes[1] with $names and " *
                "meshes[$idx] with $(keys(vertex_attributes(meshes[idx])))."
            )
        end

        consistent_face_views = true
        for name in names
            is_face_view = getproperty(m1, name) isa FaceView
            for i in 2:length(meshes)
                if (getproperty(meshes[i], name) isa FaceView) != is_face_view
                    consistent_face_views = false
                    @goto DOUBLE_BREAK
                end
            end
        end
        @label DOUBLE_BREAK

        if consistent_face_views

            # All the same kind of face, can just merge
            new_attribs = NamedTuple{names}(map(names) do name
                if name === :position
                    D = maximum(ndims, meshes)
                    T = mapreduce(m -> eltype(pointtype(m)), promote_type, meshes)
                    PT = Point{D, T}
                    return reduce(vcat, decompose.(PT, meshes))
                else
                    return reduce(vcat, getproperty.(meshes, name))
                end
            end)
            fs = reduce(vcat, faces.(meshes))

            # TODO: is the type difference in offset bad?
            idx = length(faces(m1))
            offset = length(coordinates(m1))
            views = isempty(m1.views) ? UnitRange{UInt32}[1:idx] : copy(m1.views)

            for mesh in Iterators.drop(meshes, 1)
                N = length(faces(mesh))

                # update face indices
                for i = idx .+ (1:N)
                    # TODO: face + Int changes type to Int
                    fs[i] = typeof(fs[i])(fs[i] .+ offset)
                end

                # add views
                if isempty(mesh.views)
                    push!(views, idx+1 : idx+N)
                else
                    for view in mesh.views
                        push!(views, view .+ idx)
                    end
                end

                idx += N
                offset += length(coordinates(mesh))
            end

            return Mesh(new_attribs, fs, views)

        else # mixed FaceViews and Arrays

            # simplify to VertexFace types, then retry merge
            return merge(expand_faceviews.(meshes))

        end

    end
end

"""
    expand_faceviews(mesh::Mesh)

Returns the given `mesh` if it contains no FaceViews. Otherwise, generates a new
mesh that contains no FaceViews, reordering and duplicating vertex attributes
as necessary. If the mesh has `views` they will be adjusted as needed to produce
the same submeshes.
"""
function expand_faceviews(mesh::Mesh)
    main_fs = faces(mesh)
    va = vertex_attributes(mesh)

    names = filter(name -> va[name] isa FaceView, keys(va))
    isempty(names) && return mesh

    other_fs = faces.(getproperty.((mesh,), names))
    names = (:position, names...)
    all_fs = tuple(main_fs, other_fs...)

    if isempty(mesh.views)

        new_fs, maps = merge_vertex_indices(all_fs)

        named_maps = NamedTuple{tuple(names...)}(maps)

        new_va = NamedTuple{keys(va)}(map(keys(va)) do name
            values(va[name])[get(named_maps, name, maps[1])]
        end)

        return Mesh(new_va, new_fs)

    else

        new_fs = sizehint!(eltype(main_fs)[], length(main_fs))
        new_views = sizehint!(UnitRange{UInt32}[], length(mesh.views))
        new_va = NamedTuple{keys(va)}(map(keys(va)) do name
            sizehint!(similar(values(va[name]), 0), length(va[name]))
        end)

        vertex_index_counter = eltype(first(main_fs))(1)

        for idxs in mesh.views
            view_fs, maps = merge_vertex_indices(view.(all_fs, (idxs,)), vertex_index_counter)

            vertex_index_counter += length(maps[1])

            for name in keys(new_va)
                map = maps[something(findfirst(==(name), names), 1)]
                append!(new_va[name], values(va[name])[map])
            end

            # add new faces and new view
            start = length(new_fs) + 1
            append!(new_fs, view_fs)
            push!(new_views, start:length(new_fs))
        end

        return Mesh(new_va, new_fs, new_views)
    end
end

function merge_vertex_indices(
        faces::NTuple{N_Attrib, <: AbstractVector{FT}},
        vertex_index_counter::Integer = T(1)
    ) where {N, T, FT <: AbstractFace{N, T}, N_Attrib}

    N_faces = length(faces[1])

    # maps a combination of old indices in MultiFace to a new vertex_index
    vertex_index_map = Dict{NTuple{N_Attrib, T}, T}()

    # Faces after conversion
    new_faces = sizehint!(FT[], N_faces)

    # indices that remap attributes
    attribute_indices = ntuple(n -> sizehint!(UInt32[], N_faces), N_Attrib)

    # keep track of the remapped indices for one vertex so we don't have to
    # query the dict twice
    temp = Vector{T}(undef, N)

    for multi_face in zip(faces...)

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


"""
    split_mesh(mesh::Mesh, views::Vector{UnitRange{Int}} = mesh.views)

Creates a new mesh containing `faces(mesh)[range]` for each range in `views`.
This also removes unused vertices.
"""
function split_mesh(mesh::Mesh, views::Vector{<: UnitRange{<: Integer}} = mesh.views)
    return map(views) do idxs
        new_fs, maps = merge_vertex_indices((view(faces(mesh), idxs),))

        names = keys(vertex_attributes(mesh))
        new_va = NamedTuple{names}(map(names) do name
            v = getproperty(mesh, name)
            if v isa FaceView
                _fs, _maps = merge_vertex_indices((view(faces(v), idxs),))
                return FaceView(values(v)[_maps[1]], _fs)
            else
                return v[maps[1]]
            end
        end)

        return Mesh(new_va, new_fs)
    end
end

"""
    remove_duplicates(faces)

Uses a Dict to remove duplicates from the given `faces`.
"""
function remove_duplicates(fs::AbstractVector{FT}) where {FT <: AbstractFace}
    hashmap = Dict{UInt64, FT}()
    foreach(f -> hashmap[cyclic_hash(f)] = f, fs)
    return collect(values(hashmap))
end


function map_coordinates(f, mesh::Mesh)
    result = deepcopy(mesh)
    map_coordinates!(f, result)
    return result
end

function map_coordinates(f, mesh::MetaMesh)
    result = deepcopy(Mesh(mesh))
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
        println(io, "    vertex $(name): ", attrib isa FaceView ? length(attrib.data) : length(attrib))
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
