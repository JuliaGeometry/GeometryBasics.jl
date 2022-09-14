module GeometryBasicsCore

using Base: @propagate_inbounds
using LinearAlgebra, IterTools, StaticArrays

include("offsetintegers.jl")
include("fixed_arrays.jl")
include("basic_types.jl")

include("primitives/rectangles.jl")
include("primitives/spheres.jl")
include("primitives/cylinders.jl")
include("primitives/pyramids.jl")
include("primitives/particles.jl")

include("metadata.jl")
include("interfaces.jl")
include("viewtypes.jl")
include("geometry_primitives.jl")
include("boundingboxes.jl")
include("lines.jl")

include("deprecated.jl")

# types
export AbstractGeometry, AbstractPoint, AbstractFace, GeometryPrimitive
export Mat, Point, Vec
export TupleView, FaceView, SimpleFaceView
export LineFace, Polytope, Line, NgonFace, Triangle, TriangleP
export LineString, AbstractPolygon, Polygon, MultiPoint, MultiLineString, MultiPolygon
export Simplex, Triangle, NSimplex, Tetrahedron, SimplexFace, Mesh
export TriangleFace, QuadFace, TetrahedronFace
export Tesselation, Normal, UV, UVW
export HyperRectangle, Rect, Rect2, Rect3, Recti, Rect2i, Rect3i, Rectf, Rect2f, Rect3f
export HyperSphere, Circle, Sphere, Cylinder, Cylinder2, Cylinder3, Pyramid
export OffsetInteger, ZeroIndex, OneIndex, GLIndex
export MetaT

# methods
export max_dist_dim, max_euclidean, max_euclideansq, min_dist_dim, min_euclidean
export min_euclideansq, minmax_dist_dim, minmax_euclidean, minmax_euclideansq
export connect, extremity, coordinates, convert_simplex
export height, origin, radius, width, widths
export centered, direction, area, volume, update
export faces, normals, texturecoordinates
export before, during, meets, overlaps, intersects, finishes
export self_intersections, split_intersections
export decompose, decompose_uv, decompose_normals
export meta, metafree

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end
