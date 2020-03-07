const FaceMesh{Dim, T, Element} = Mesh{Dim, T, Element, <: FaceView{Element}}
coordinates(mesh::FaceMesh) = coordinates(getfield(mesh, :simplices))
faces(mesh::FaceMesh) = faces(getfield(mesh, :simplices))

const GLTriangleElement = Triangle{3, Float32}
const GLTriangleFace = TriangleFace{GLIndex}

const TriangleMesh{Dim, T, PointType} = AbstractVector{TriangleP{Dim, T, PointType}}

const PlainTriangleMesh{Dim, T} = TriangleMesh{Dim, T, Point{Dim, T}}

const GLPlainTriangleMesh = PlainTriangleMesh{3, Float32}

const PointWithUV{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:uv,), Tuple{Vec{2, T}}}
const PointWithNormal{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:normals,), Tuple{Vec{3, T}}}
const PointWithUVNormal{Dim, T} = PointMeta{Dim, T, Point{Dim, T}, (:normals, :uv), Tuple{Vec{3, T}, Vec{2, T}}}

const UVMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUV{Dim, T}}
const GLUVMesh{Dim} = UVMesh{Dim, Float32}
const GLUVMesh2D = UVMesh{2}
const GLUVMesh3D = UVMesh{3}

const UVNormalMesh{Dim, T} = TriangleMesh{Dim, T, PointWithUVNormal{Dim, T}}
const GLUVNormalMesh{Dim} = UVNormalMesh{Dim, Float32}
const GLUVNormalMesh2D = UVNormalMesh{2}
const GLUVNormalMesh3D = UVNormalMesh{3}

const NormalMesh{Dim, T} = TriangleMesh{Dim, T, PointWithNormal{Dim, T}}
const GLNormalMesh{Dim} = NormalMesh{Dim, Float32}
const GLNormalMesh2D = GLNormalMesh{2}
const GLNormalMesh3D = GLNormalMesh{3}


function triangle_mesh(primitive)
    return Mesh(decompose(Point3f0, primitive), decompose(GLTriangleFace, primitive))
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
    return Mesh(meta(points; uv=uv, normals=normals(points, fs)), fs)
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
