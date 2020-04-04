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
Abstract Mesh with triangle elements
"""
const TriangleMesh{Dim, T, PointType} = AbstractMesh{TriangleP{Dim, T, PointType}}

"""
    PlainMesh
Triangle mesh with no meta information (just points + triangle faces)
"""
const PlainMesh{Dim, T} = TriangleMesh{Dim, T, Point{Dim, T}}
const GLPlainMesh{Dim} = PlainMesh{Dim, Float32}
const GLPlainMesh2D = GLPlainMesh{2}
const GLPlainMesh3D = GLPlainMesh{3}

"""
    UVMesh
PlainMesh with texture coordinates meta at each point.
`uvmesh.uv isa AbstractVector{Vec2f0}`
"""
const UVMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUV{Dim, T}}
const GLUVMesh{Dim} = UVMesh{Dim, Float32}
const GLUVMesh2D = UVMesh{2}
const GLUVMesh3D = UVMesh{3}

"""
    NormalMesh
PlainMesh with normals meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
"""
const NormalMesh{Dim, T} = TriangleMesh{Dim, T, PointWithNormal{Dim, T}}
const GLNormalMesh{Dim} = NormalMesh{Dim, Float32}
const GLNormalMesh2D = GLNormalMesh{2}
const GLNormalMesh3D = GLNormalMesh{3}

"""
    NormalUVMesh
PlainMesh with normals and uv meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
`normalmesh.uv isa AbstractVector{Vec2f0}`
"""
const NormalUVMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUVNormal{Dim, T}}
const GLNormalUVMesh{Dim} = NormalUVMesh{Dim, Float32}
const GLNormalUVMesh2D = GLNormalUVMesh{2}
const GLNormalUVMesh3D = GLNormalUVMesh{3}

"""
    NormalUVWMesh
PlainMesh with normals and uvw (texture coordinates in 3D) meta at each point.
`normalmesh.normals isa AbstractVector{Vec3f0}`
`normalmesh.uvw isa AbstractVector{Vec3f0}`
"""
const NormalUVWMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUVWNormal{Dim, T}}
const GLNormalUVWMesh{Dim} = NormalUVWMesh{Dim, Float32}
const GLNormalUVWMesh2D = GLNormalUVWMesh{2}
const GLNormalUVWMesh3D = GLNormalUVWMesh{3}

function triangle_mesh(primitive, nvertices=25)
    return Mesh(decompose(Point3f0, primitive, nvertices), decompose(GLTriangleFace, primitive, nvertices))
end

function triangle_mesh(primitive::AbstractVector{<:AbstractPoint})
    return Mesh(primitive, decompose(GLTriangleFace, primitive))
end

function triangle_mesh(primitive::GeometryPrimitive{N, T}) where {N, T}
    return Mesh(decompose(Point{N, T}, primitive), decompose(GLTriangleFace, primitive))
end

function gl_triangle_mesh(primitive::GeometryPrimitive{N}) where {N}
    return Mesh(decompose(Point{N, Float32}, primitive), decompose(GLTriangleFace, primitive))
end

function gl_triangle_mesh(mesh::GLPlainMesh)
    return mesh
end

function gl_triangle_mesh(poly::AbstractVector{<: AbstractPoint})
    return Mesh(decompose(Point{2, Float32}, poly), decompose(GLTriangleFace, poly))
end

function uv_triangle_mesh(primitive::GeometryPrimitive{N, T}) where {N, T}
    points = decompose(Point{N, T}, primitive)
    fs = decompose(GLTriangleFace, primitive)
    uv = collect_with_eltype(Vec2f0, texturecoordinates(primitive))
    return Mesh(meta(points; uv=uv), fs)
end

function gl_uv_triangle_mesh3d(primitive::GeometryPrimitive)
    points = decompose(Point3f0, primitive)
    fs = decompose(GLTriangleFace, primitive)
    uv = decompose_uv(primitive)
    return Mesh(meta(points; uv=uv), fs)
end

function gl_uvw_triangle_mesh3d(primitive::GeometryPrimitive)
    points = decompose(Point3f0, primitive)
    fs = decompose(GLTriangleFace, primitive)
    uvw = decompose_uvw(primitive)
    return Mesh(meta(points; uvw=uvw), fs)
end

function gl_uv_triangle_mesh2d(primitive::GeometryPrimitive)
    points = decompose(Point2f0, primitive)
    fs = decompose(GLTriangleFace, primitive)
    uv = decompose_uv(primitive)
    return Mesh(meta(points; uv=uv), fs)
end

function gl_normal_mesh3d(primitive::GeometryPrimitive)
    points = decompose(Point3f0, primitive)
    fs = decompose(GLTriangleFace, primitive)
    return Mesh(meta(points; normals=normals(points, fs)), fs)
end

function normal_mesh(primitive::GeometryPrimitive{N, T}) where {N, T}
    points = decompose(Point{N, T}, primitive)
    fs = decompose(GLTriangleFace, primitive)
    return Mesh(meta(points; normals=normals(points, fs)), fs)
end

"""
Calculate the signed volume of one tetrahedron. Be sure the orientation of your surface is right.
"""
function volume(triangle::Triangle) where {VT,FT}
    v1, v2, v3 = triangle
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ ( v2 × v3 )) / 6
end

"""
Calculate the signed volume of all tetrahedra. Be sure the orientation of your surface is right.
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

function decompose_normals(mesh::Mesh)
    normal_vectors = normals(mesh)
    return convert(Vector{Vec3f0}, normal_vectors)
end
