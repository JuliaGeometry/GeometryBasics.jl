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

    if positions isa FaceView
        positions = positions.data
        _fs = positions.faces
    else
        _fs = faces(primitive)
    end

    # If faces returns nothing for primitive, we try to triangulate!
    if isnothing(_fs)
        if eltype(positions) <: Point2
            # triangulation.jl
            fs = decompose(facetype, positions)
        else
            error("No triangulation for $(typeof(primitive))")
        end
    else
        fs = _fs
    end

    return mesh(positions, collect(fs); facetype = facetype, vertex_attributes...)
end

function drop_nothing_kwargs(kwargs)
    _keys = filter(k -> !isnothing(kwargs[k]), keys(kwargs))
    _vals = ntuple(n -> kwargs[_keys[n]], length(_keys))
    return NamedTuple{_keys}(_vals)
end

"""
    mesh(positions, faces[, facetype = GLTriangleFace, vertex_attributes...])

Creates a mesh from the given positions and faces. Other vertex attributes like
normals and texture coordinates can be added as keyword arguments.

By default the generated mesh uses `GLTriangleFace`. If the input faces are of 
type `MultiFace` they will get converted appropriately, which may cause 
reordering and duplication of positions and other vertex attributes.
"""
function mesh(
        positions::AbstractVector{<:Point}, 
        faces::AbstractVector{FT}; 
        facetype=GLTriangleFace, vertex_attributes...
    ) where {FT <: AbstractVertexFace}

    fs = decompose(facetype, faces)
    va = drop_nothing_kwargs(vertex_attributes)

    if (facetype != FT)
        va = NamedTuple(map(keys(va)) do name
            attrib = va[name]
            if attrib isa FaceView
                return name => FaceView(attrib.data, decompose(facetype, attrib.faces))
            else
                return name => attrib
            end
        end)
    end

    return Mesh(positions, fs; va...)
end

function mesh(
        mesh::Mesh{D, T, FT}; pointtype = Point{D, Float32}, 
        facetype::Type{<: AbstractFace} = GLTriangleFace,
        attributes...
    ) where {D, T, FT <: AbstractFace}
    
    va = drop_nothing_kwargs(attributes)

    # N = length(mesh.position)
    # if !all(attr -> length(attr) == N, values(va))
    #     error("At least one of the given vertex attributes does not match `length(mesh.positon) = $N`.")
    # end

    if FT == facetype
        if isempty(va) && GeometryBasics.pointtype(mesh) == pointtype
            return mesh
        else
            # 1. add vertex attributes, 2. convert position attribute
            va = merge(vertex_attributes(mesh), va)
            va = merge(va, (position = decompose(pointtype, va.position),))
            return Mesh(va, faces(mesh), mesh.views)
        end
    else
        # 1. add vertex attributes, 2. convert position attribute
        va = merge(vertex_attributes(mesh), va)
        va = merge(va, (position = decompose(pointtype, va.position),))

        # Resolve facetype conversions of FaceViews
        va = NamedTuple(map(keys(va)) do name
            attrib = va[name]
            if attrib isa FaceView
                return name => FaceView(attrib.data, decompose(facetype, attrib.faces))
            else
                return name => attrib
            end
        end)

        # update main face type
        f, views = decompose(facetype, faces(mesh), mesh.views)
        return Mesh(va, f, views)
    end
end

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)

Create a mesh from a polygon given as a vector of points, using triangulation.
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace) where {P<:Point{2}}
    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype)
end

function triangle_mesh(
        primitive::Union{AbstractGeometry{N}, AbstractVector{<: Point{N}}};
        pointtype = Point{N, Float32}, facetype = GLTriangleFace
    )::SimpleTriangleMesh{N} where {N} # TODO: is the output type doing anything here?

    return mesh(primitive; pointtype = pointtype, facetype = facetype)
end


pointtype(::Mesh{D, T}) where {D, T} = Point{D, T}
facetype(::Mesh{D, T, FT}) where {D, T, FT} = FT

function triangle_mesh(primitive::Mesh{N}; pointtype = Point{N, Float32}) where {N}
    # already target type:
    if GeometryBasics.pointtype(primitive) === pointtype && 
            GeometryBasics.facetype(primitive) === GLTriangleFace
        return primitive
    else
        return mesh(primitive; pointtype = pointtype, facetype = GLTriangleFace)
    end
end


function uv_mesh(
        primitive::AbstractGeometry{N}; pointtype = Point{N, Float32}, 
        uvtype = Vec2f, facetype = GLTriangleFace
    ) where {N}
    
    return mesh(
        primitive, uv = decompose_uv(uvtype, primitive), 
        pointtype = pointtype, facetype = facetype
    )
end

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

function normal_mesh(
        points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractVertexFace}; 
        pointtype = Point3f, normaltype = Vec3f, facetype = GLTriangleFace
    )
    _points = decompose(pointtype, points)
    _faces = decompose(facetype, faces)
    return Mesh(_faces, position = _points, normal = normals(_points, _faces, normaltype))
end

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
                return mapreduce(m -> getproperty(m, name), vcat, meshes)
            end)
            fs = reduce(vcat, faces.(meshes))

            # TODO: is the type difference in offset bad?
            idx = length(faces(m1))
            offset = length(coordinates(m1))
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
                offset += length(coordinates(mesh))
            end

            return Mesh(new_attribs, fs, views)

        else # mixed FaceViews and Arrays

            # simplify to VertexFace types, then retry merge
            return merge(clear_faceviews.(meshes))

        end

    end
end

function clear_faceviews(mesh::Mesh)
    main_fs = faces(mesh)
    va = vertex_attributes(mesh)

    names = filter(name -> va[name] isa FaceView, collect(keys(va)))
    isempty(names) && return mesh

    other_fs = map(name -> va[name].faces, names)
    pushfirst!(names, :position)
    all_fs = tuple(main_fs, other_fs...)
    
    
    if isempty(mesh.views)

        new_fs, maps = merge_vertex_indices(all_fs)

        named_maps = NamedTuple{tuple(names...)}(maps)

        new_va = NamedTuple(map(collect(keys(va))) do name
            attrib = va[name]
            if name === :position
                return name => attrib[named_maps[name]]
            elseif haskey(named_maps, name)
                return name => attrib.data[named_maps[name]]
            else
                return name => attrib
            end
        end)

        return Mesh(new_va, new_fs)

    else

        new_fs = sizehint!(eltype(main_fs)[], length(main_fs))
        new_views = sizehint!(UnitRange{Int}[], length(mesh.views))
        new_va = NamedTuple(map(collect(keys(va))) do name
            attrib = va[name]
            if name == :position
                return name => sizehint!(eltype(attrib)[], length(attrib))
            elseif name in names
                return name => sizehint!(eltype(attrib.data)[], length(attrib.data))
            else # doesn't need remapping so we just add it here
                return name => attrib
            end
        end)

        vertex_index_counter = eltype(first(main_fs))(1)

        for idxs in mesh.views
            view_fs, maps = merge_vertex_indices(view.(all_fs, (idxs,)), vertex_index_counter)

            vertex_index_counter += length(maps[1])
            
            for (name, map) in zip(names, maps)
                if name == :position
                    append!(new_va[name], va[name][map])
                else
                    append!(new_va[name], va[name].data[map])
                end
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
        faces::AbstractVector{FT}, args...
    ) where {N, T, FT <: AbstractFace{N, T}}
    if args[end] isa Integer
        fs = tuple(faces, args[1:end-1]...)
        return merge_vertex_indices(fs, args[end])
    else
        return merge_vertex_indices(tuple(faces, args...))
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

    # keep track of the remmaped indices for one vertex so we don't have to 
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
