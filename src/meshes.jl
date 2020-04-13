const FaceMesh{Dim, T, Element} = Mesh{Dim, T, Element, <: FaceView{Element}}

coordinates(mesh::FaceMesh) = coordinates(getfield(mesh, :simplices))
faces(mesh::FaceMesh) = faces(getfield(mesh, :simplices))

const GLTriangleElement = Triangle{3, Float32}
const GLTriangleFace = TriangleFace{GLIndex}
const PointWithUV{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:uv,), Tuple{Vec{2, T}}}
const PointWithNormal{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:normals,), Tuple{Vec{3, T}}}
const PointWithUVNormal{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:normals, :uv), Tuple{Vec{3, T}, Vec{2, T}}}
const PointWithUVWNormal{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:normals, :uvw), Tuple{Vec{3, T}, Vec{3, T}}}

"""
    TriangleMesh{Dim, T, PointType}

Abstract Mesh with triangle elements of eltype `T`.
"""
const TriangleMesh{Dim, T, PointType} = AbstractMesh{TriangleP{Dim, T, PointType}}

"""
    PlainMesh{Dim, T}

Triangle mesh with no meta information (just points + triangle faces)
"""
const PlainMesh{Dim, T} = TriangleMesh{Dim, T, Point{Dim, T}}
const GLPlainMesh{Dim} = PlainMesh{Dim, Float32}
const GLPlainMesh2D = GLPlainMesh{2}
const GLPlainMesh3D = GLPlainMesh{3}

"""
    UVMesh{Dim, T}

PlainMesh with texture coordinates meta at each point.
`uvmesh.uv isa AbstractVector{Vec2f0}`
"""
const UVMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUV{Dim, T}}
const GLUVMesh{Dim} = UVMesh{Dim, Float32}
const GLUVMesh2D = UVMesh{2}
const GLUVMesh3D = UVMesh{3}

"""
    NormalMesh{Dim, T}

PlainMesh with normals meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
"""
const NormalMesh{Dim, T} = TriangleMesh{Dim, T, PointWithNormal{Dim, T}}
const GLNormalMesh{Dim} = NormalMesh{Dim, Float32}
const GLNormalMesh2D = GLNormalMesh{2}
const GLNormalMesh3D = GLNormalMesh{3}

"""
    NormalUVMesh{Dim, T}

PlainMesh with normals and uv meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
`normalmesh.uv isa AbstractVector{Vec2f0}`
"""
const NormalUVMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUVNormal{Dim, T}}
const GLNormalUVMesh{Dim} = NormalUVMesh{Dim, Float32}
const GLNormalUVMesh2D = GLNormalUVMesh{2}
const GLNormalUVMesh3D = GLNormalUVMesh{3}

"""
    NormalUVWMesh{Dim, T}

PlainMesh with normals and uvw (texture coordinates in 3D) meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
`normalmesh.uvw isa AbstractVector{Vec3f0}`
"""
const NormalUVWMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUVWNormal{Dim, T}}
const GLNormalUVWMesh{Dim} = NormalUVWMesh{Dim, Float32}
const GLNormalUVWMesh2D = GLNormalUVWMesh{2}
const GLNormalUVWMesh3D = GLNormalUVWMesh{3}

best_nvertices(any_primitive) = 24
best_pointtype(::GeometryPrimitive{N, T}) where {N, T} = Point{N, T}


"""
    mesh(primitive::GeometryPrimitive;
         pointtype=best_pointtype(primitive), facetype=GLTriangle,
         uvtype=nothing, normaltype=nothing, nvertices=nothing)

Creates a mesh from `primitive`.
Uses the element types from the keyword arguments to create the attributes.
The attributes that have their type set to nothing are not added to the mesh.
`nvertices` specifies the numbers of vertices the primitive is decomposed to.
Note, that this can be an `Int` or `Tuple{Int, Int}``, when the primitive is grid based.
It also only losely correlates to the number of vertices, depending on the algorithm used.
#TODO: find a better number here!
"""
function mesh(primitive::GeometryPrimitive;
              pointtype=best_pointtype(primitive), facetype=GLTriangleFace,
              uv=nothing, normaltype=nothing, nvertices=nothing)

    if nvertices === nothing
        nvertices = best_nvertices(primitive)
    end

    positions = decompose(pointtype, primitive, nvertices)
    faces = decompose(facetype, primitive, nvertices)
    # If faces returns nothing for primitive, we try to triangulate!
    if faces === nothing
        # triangulation.jl
        faces = decompose(facetype, positions)
    end
    attributes = Dict{Symbol, Any}()

    if uv !== nothing
        attributes[:uv] = decompose(UV(uv), primitive, nvertices)
    end

    if normaltype !== nothing
        primitive_normals = decompose(Normal(normaltype), primitive, nvertices)
        if primitive_normals !== nothing
            attributes[:normals] = primitive_normals
        else
            # Normals not implemented for primitive, so we calculate them!
            attributes[:normals] = normals(positions, faces)
        end
    end

    return Mesh(meta(positions; attributes...), faces)
end

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)
Polygon triangluation!
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
              normaltype=nothing, nvertices=nothing) where {P<:AbstractPoint{2}}

    faces = decompose(facetype, polygon)
    positions = decompose(pointtype, polygon)

    if nvertices !== nothing
        error("Resampling polygon not supported!")
    end

    if normaltype !== nothing
        n = normals(positions, faces; normaltype=normaltype)
        positions = meta(positions; normals=n)
    end

    return Mesh(positions, faces)
end

mesh(m::Mesh) = m
triangle_mesh(m::Mesh) = m

function triangle_mesh(primitive::GeometryPrimitive{N}; nvertices=nothing) where {N}
    return mesh(primitive; pointtype=Point{N, Float32}, facetype=GLTriangleFace, nvertices=nvertices)
end


function triangle_mesh(polygon::AbstractVector{<: AbstractPoint{2}})
    return mesh(polygon; pointtype=Point{2, Float32}, facetype=GLTriangleFace)
end

function uv_mesh(primitive::GeometryPrimitive{N, T}; nvertices=nothing) where {N, T}
    return mesh(primitive; pointtype=Point{N, Float32}, uv=Vec2f0,
                facetype=GLTriangleFace, nvertices=nvertices)
end

function uv_normal_mesh(primitive::GeometryPrimitive{N}; nvertices=nothing) where {N}
    return mesh(primitive; pointtype=Point{N, Float32}, uv=Vec2f0, normaltype=Vec3f0,
                facetype=GLTriangleFace, nvertices=nvertices)
end

function normal_mesh(primitive::GeometryPrimitive{N}; nvertices=nothing) where {N}
    return mesh(primitive; pointtype=Point{N, Float32}, normaltype=Vec3f0,
                facetype=GLTriangleFace, nvertices=nvertices)
end

function normal_mesh(points::AbstractVector{<:AbstractPoint},
                     faces::AbstractVector{<:AbstractFace})
    _points = convert(Vector{Point3f0}, points)
    _faces = decompose(GLTriangleFace, faces)
    return Mesh(meta(_points; normals=normals(_points, _faces)), _faces)
end

"""
    volume(triangle)

Calculate the signed volume of one tetrahedron. Be sure the orientation of your
surface is right.
"""
function volume(triangle::Triangle) where {VT,FT}
    v1, v2, v3 = triangle
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ ( v2 × v3 )) / 6
end

"""
    volume(mesh)

Calculate the signed volume of all tetrahedra. Be sure the orientation of your
surface is right.
"""
function volume(mesh::Mesh) where {VT, FT}
    return sum(volume, mesh)
end

function Base.merge(meshes::AbstractVector{<: Mesh})
    if isempty(meshes)
        error("No meshes to merge")
    elseif length(meshes) == 1
        return meshes[1]
    else
        m1 = meshes[1]
        ps = copy(coordinates(m1))
        fs = copy(faces(m1))
        for mesh in Iterators.drop(meshes, 1)
            append!(fs, map(f-> f .+ length(ps), faces(mesh)))
            append!(ps, coordinates(mesh))
        end
        return Mesh(ps, fs)
    end
end

function decompose(::Normal{T}, mesh::Mesh) where {T}
    normal_vectors = normals(mesh)
    return decompose(T, normal_vectors)
end

function decompose(::UV{T}, mesh::Mesh) where {T}
    if hasproperty(mesh, :uv)
        return decompose(T, mesh.uv)
    else
        error("Mesh doesn't have UV texture coordinates")
    end
end

"""
    pointmeta(mesh::Mesh; meta_data...)

Attaches metadata to the coordinates of a mesh
"""
function pointmeta(mesh::Mesh; meta_data...)
    points = coordinates(mesh)
    attr = GeometryBasics.attributes(points)
    delete!(attr, :position) # position == metafree(points)
    # delete overlapping attributes so we can replace with `meta_data`
    foreach(k-> delete!(attr, k), keys(meta_data))
    return Mesh(meta(metafree(points); attr..., meta_data...), faces(mesh))
end

"""
    pop_pointmeta(mesh::Mesh, property::Symbol)
Remove `property` from point metadata.
Returns the new mesh, and the property!
"""
function pop_pointmeta(mesh::Mesh, property::Symbol)
    points = coordinates(mesh)
    attr = GeometryBasics.attributes(points)
    delete!(attr, :position) # position == metafree(points)
    # delete overlapping attributes so we can replace with `meta_data`
    m = pop!(attr, property)
    return Mesh(meta(metafree(points); attr...), faces(mesh)), m
end

"""
    facemeta(mesh::Mesh; meta_data...)

Attaches metadata to the faces of a mesh
"""
function facemeta(mesh::Mesh; meta_data...)
    return Mesh(coordinates(mesh), meta(faces(mesh); meta_data...))
end
