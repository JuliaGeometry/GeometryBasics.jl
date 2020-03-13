module GeometryBasics

    using StaticArrays, Tables, StructArrays, IterTools, LinearAlgebra

    using Base: @propagate_inbounds

    include("fixed_arrays.jl")
    include("offsetintegers.jl")
    include("basic_types.jl")
    include("metadata.jl")
    include("viewtypes.jl")
    include("geometry_primitives.jl")
    include("rectangles.jl")
    include("triangulation.jl")
    include("meshes.jl")
    include("lines.jl")
    include("boundingboxes.jl")

    export AbstractGeometry, GeometryPrimitive
    export Mat, Point, Vec
    export LineFace, Polytope, Line, NgonFace, simplex_convert
    export LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon
    export Simplex, connect, Triangle, NSimplex, Tetrahedron
    export QuadFace, metafree, coordinates, TetrahedronFace
    export TupleView, SimplexFace, Mesh, meta
    export Triangle, TriangleP
    export AbstractFace, TriangleFace, QuadFace, GLTriangleFace
    export OffsetInteger, ZeroIndex, OneIndex, GLIndex
    export FaceView, SimpleFaceView
    export AbstractPoint, PointMeta, PointWithUV
    export decompose, coordinates, faces, normals
    export GLTriangleFace, GLNormalMesh3D, GLPlainTriangleMesh, GLUVMesh3D, GLUVNormalMesh3D
    export AbstractMesh, Mesh, TriangleMesh
    export GLNormalMesh2D, PlainTriangleMesh

    # all the different predefined mesh types
    # Note: meshes can contain arbitrary meta information,
    export AbstractMesh, TriangleMesh, PlainMesh, GLPlainMesh, GLPlainMesh2D, GLPlainMesh3D
    export UVMesh, GLUVMesh, GLUVMesh2D, GLUVMesh3D
    export NormalMesh, GLNormalMesh, GLNormalMesh2D, GLNormalMesh3D
    export NormalUVMesh, GLNormalUVMesh, GLNormalUVMesh2D, GLNormalUVMesh3D
    export NormalUVWMesh, GLNormalUVWMesh, GLNormalUVWMesh2D, GLNormalUVWMesh3D

    # mesh creation functions
    export triangle_mesh, gl_triangle_mesh, uv_triangle_mesh
    export gl_uv_triangle_mesh3d, gl_normal_mesh3d, normal_mesh

    export height, origin, radius, width, widths, xwidth, yheight
    export HyperSphere, Circle, Sphere
    export Cylinder, Cylinder2, Cylinder3, Pyramid, extremity
    export Rect, Rect2D, Rect3D, IRect, IRect2D, IRect3D, FRect, FRect2D, FRect3D
    export before, during, contains, isinside, isoutside, meets, overlaps, intersects, finishes
    export centered, direction, area, update
    export max_dist_dim, max_euclidean, max_euclideansq, min_dist_dim, min_euclidean
    export min_euclideansq, minmax_dist_dim, minmax_euclidean, minmax_euclideansq
    export self_intersections, split_intersections

end # module
