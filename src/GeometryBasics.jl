module GeometryBasics

using StaticArrays, Tables, StructArrays, IterTools, LinearAlgebra
using GeoInterface
import Extents
using EarCut_jll

using Base: @propagate_inbounds

include("fixed_arrays.jl")
include("offsetintegers.jl")
include("basic_types.jl")

include("primitives/rectangles.jl")
include("primitives/spheres.jl")
include("primitives/cylinders.jl")
include("primitives/pyramids.jl")
include("primitives/particles.jl")

include("interfaces.jl")
include("metadata.jl")
include("viewtypes.jl")
include("geometry_primitives.jl")
include("meshes.jl")
include("triangulation.jl")
include("lines.jl")
include("boundingboxes.jl")

include("deprecated.jl")
include("geointerface.jl")

export AbstractGeometry, GeometryPrimitive
export Mat, Point, Vec
export LineFace, Polytope, Line, NgonFace, convert_simplex
export LineString, AbstractPolygon, Polygon, MultiPoint, MultiLineString, MultiPolygon
export Simplex, connect, Triangle, NSimplex, Tetrahedron
export QuadFace, metafree, coordinates, TetrahedronFace
export TupleView, SimplexFace, Mesh, meta
export Triangle, TriangleP
export AbstractFace, TriangleFace, QuadFace, GLTriangleFace
export OffsetInteger, ZeroIndex, OneIndex, GLIndex
export FaceView, SimpleFaceView
export AbstractPoint, PointMeta, PointWithUV
export PolygonMeta, MultiPointMeta, MultiLineStringMeta, MeshMeta, LineStringMeta,
       MultiPolygonMeta
export decompose, coordinates, faces, normals, decompose_uv, decompose_normals,
       texturecoordinates
export Tesselation, pointmeta, Normal, UV, UVW
export GLTriangleFace, GLUVMesh3D
export AbstractMesh, Mesh, TriangleMesh
export GLNormalMesh2D
export MetaT, meta_table

# all the different predefined mesh types
# Note: meshes can contain arbitrary meta information,
export AbstractMesh, TriangleMesh, PlainMesh, GLPlainMesh, GLPlainMesh2D, GLPlainMesh3D
export UVMesh, GLUVMesh, GLUVMesh2D, GLUVMesh3D
export NormalMesh, GLNormalMesh, GLNormalMesh2D, GLNormalMesh3D
export NormalUVMesh, GLNormalUVMesh, GLNormalUVMesh2D, GLNormalUVMesh3D
export NormalUVWMesh, GLNormalUVWMesh, GLNormalUVWMesh2D, GLNormalUVWMesh3D

# mesh creation functions
export triangle_mesh, triangle_mesh, uv_mesh
export uv_mesh, normal_mesh, uv_normal_mesh

export height, origin, radius, width, widths
export HyperSphere, Circle, Sphere
export Cylinder, Cylinder2, Cylinder3, Pyramid, extremity
export HyperRectangle, Rect, Rect2, Rect3, Recti, Rect2i, Rect3i, Rectf, Rect2f, Rect3f
export before, during, meets, overlaps, intersects, finishes
export centered, direction, area, volume, update
export max_dist_dim, max_euclidean, max_euclideansq, min_dist_dim, min_euclidean
export min_euclideansq, minmax_dist_dim, minmax_euclidean, minmax_euclideansq
export self_intersections, split_intersections

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end # module
