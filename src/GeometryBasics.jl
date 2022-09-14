module GeometryBasics

using Base: @propagate_inbounds

using Reexport
@reexport using GeometryBasicsCore
import GeometryBasicsCore: AbstractMesh, AbstractSimplexFace, AbstractNgonFace
import GeometryBasicsCore: VecTypes, Meshable, FaceMesh, Ngon, MetaType, MetaFree

using Tables, StructArrays, LinearAlgebra
using GeoInterface
using EarCut_jll

include("metadata.jl")
include("meshes.jl")
include("triangulation.jl")
include("geointerface.jl")

Tables.schema(mesh::Mesh) = Tables.schema(getfield(mesh, :simplices))
Tables.schema(faceview::FaceView) = Tables.schema(getfield(faceview, :elements))

@inline function GeometryBasicsCore.connect(points::AbstractMatrix{T},
                         P::Type{<:AbstractPoint{N}}) where {T <: Real,N}
    return if size(points, 1) === N
        return reinterpret(Point{N,T}, points)
    elseif size(points, 2) === N
        seglen = size(points, 1)
        columns = ntuple(N) do i
            return view(points, ((i - 1) * seglen + 1):(i * seglen))
        end
        return StructArray{Point{N,T}}(columns)
    else
        error("Dim 1 or 2 must be equal to the point dimension!")
    end
end

# types
export PointMeta, PointWithUV
export PolygonMeta, MultiPointMeta, MultiLineStringMeta, MeshMeta, LineStringMeta, MultiPolygonMeta
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
