module GeometryBasics

using Base: @propagate_inbounds

using Reexport
@reexport using GeometryBasicsCore
import GeometryBasicsCore: AbstractMesh, AbstractSimplexFace, AbstractNgonFace
import GeometryBasicsCore: VecTypes, Meshable, FaceMesh, Ngon, MetaType, MetaFree
import GeometryBasicsCore: orthogonal_vector, getnamestypes

using Tables, StructArrays, LinearAlgebra
using GeoInterface, EarCut_jll

include("core.jl")
include("metadata.jl")
include("meshes.jl")
include("triangulation.jl")
include("geointerface.jl")

# types
export PointMeta, PointWithUV
export PolygonMeta, MultiPointMeta, MultiLineStringMeta
export MeshMeta, LineStringMeta, MultiPolygonMeta
export GLTriangleFace, GLUVMesh3D, GLNormalMesh2D

# all the different predefined mesh types
# Note: meshes can contain arbitrary meta information,
export TriangleMesh, PlainMesh, GLPlainMesh, GLPlainMesh2D, GLPlainMesh3D
export UVMesh, GLUVMesh, GLUVMesh2D, GLUVMesh3D
export NormalMesh, GLNormalMesh, GLNormalMesh2D, GLNormalMesh3D
export NormalUVMesh, GLNormalUVMesh, GLNormalUVMesh2D, GLNormalUVMesh3D
export NormalUVWMesh, GLNormalUVWMesh, GLNormalUVWMesh2D, GLNormalUVWMesh3D

# methods
export triangle_mesh, uv_mesh, normal_mesh, uv_normal_mesh
export meta_table, pointmeta

if Base.VERSION >= v"1.4.2"
    include("precompile.jl")
    _precompile_()
end

end # module
