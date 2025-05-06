module GeometryBasics

using IterTools, LinearAlgebra, StaticArrays
using GeoInterface
import Extents
using EarCut_jll
import Base: *

using Base: @propagate_inbounds

include("fixed_arrays.jl")
include("offsetintegers.jl")
include("basic_types.jl")

include("primitives/rectangles.jl")
include("primitives/spheres.jl")
include("primitives/cylinders.jl")
include("primitives/pyramids.jl")
include("primitives/particles.jl")
include("primitives/Cone.jl")

include("interfaces.jl")
include("viewtypes.jl")
include("geometry_primitives.jl")
include("meshes.jl")
include("triangulation.jl")
include("lines.jl")
include("boundingboxes.jl")

include("geointerface.jl")

export AbstractGeometry, GeometryPrimitive
export Mat, Point, Vec
export LineFace, Polytope, Line, NgonFace, convert_simplex
export LineString, MultiLineString, MultiPoint
export AbstractPolygon, Polygon, MultiPolygon
export Simplex, connect, Triangle, NSimplex, Tetrahedron
export QuadFace, coordinates, TetrahedronFace
export TupleView, SimplexFace
export Triangle
export AbstractFace, TriangleFace, QuadFace, GLTriangleFace
export OffsetInteger, ZeroIndex, OneIndex, GLIndex
export decompose, coordinates, faces, normals, decompose_uv, decompose_normals,
       texturecoordinates, vertex_attributes
export expand_faceviews, split_mesh, remove_duplicates
export face_normals
export Tessellation, Normal, UV, UVW
export AbstractMesh, Mesh, MetaMesh, FaceView
export per_face


# all the different predefined mesh types
# Note: meshes can contain arbitrary meta information,
export AbstractMesh

# mesh creation functions
export triangle_mesh, triangle_mesh, uv_mesh
export uv_mesh, normal_mesh, uv_normal_mesh

export height, origin, radius, width, widths
export HyperSphere, Circle, Sphere, Cone
export Cylinder, Pyramid, extremity
export HyperRectangle, Rect, Rect2, Rect3, Recti, Rect2i, Rect3i, Rectf, Rect2f, Rect3f, Rectd, Rect2d, Rect3d, RectT
export before, during, meets, overlaps, intersects, finishes
export centered, direction, area, volume, update
export max_dist_dim, max_euclidean, max_euclideansq, min_dist_dim, min_euclidean
export min_euclideansq, minmax_dist_dim, minmax_euclidean, minmax_euclideansq
export self_intersections, split_intersections

if Base.VERSION >= v"1.8"
    include("precompiles.jl")
end

end # module
