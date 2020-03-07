

const ZeroIndex{T <: Integer} = OffsetInteger{-1, T}
const OneIndex{T <: Integer} = OffsetInteger{0, T}
const GLIndex = ZeroIndex{Cuint}

const Triangle{T} = Face{3, T}
const GLFace{Dim} = Face{Dim, GLIndex} #offset is relative to julia, so -1 is 0-indexed
const GLTriangle = Face{3, GLIndex}
const GLQuad = Face{4, GLIndex}



const VecTypes{N, T} = Union{StaticVector{N, T}, NTuple{N, T}}

const HMesh = HomogenousMesh


"""
A `SimpleMesh` is an alias for a `HomogenousMesh` parameterized only by
vertex and face types.
"""
const SimpleMesh{VT, FT} = HMesh{VT, FT, Nothing, Nothing, Nothing, Nothing, Nothing}
const PlainMesh{VT, FT} = HMesh{Point{3, VT}, FT, Nothing, Nothing, Nothing, Nothing, Nothing}
const GLPlainMesh = PlainMesh{Float32, GLTriangle}

const Mesh2D{VT, FT} = HMesh{Point{2, VT}, FT, Nothing, Nothing, Nothing, Nothing, Nothing}
const GLMesh2D = Mesh2D{Float32, GLTriangle}

const UVMesh{VT, FT, UVT} = HMesh{Point{3, VT}, FT, Nothing, UV{UVT}, Nothing, Nothing, Nothing}
const GLUVMesh = UVMesh{Float32, GLTriangle, Float32}

const UVWMesh{VT, FT, UVT} = HMesh{Point{3, VT}, FT, Nothing, UVW{UVT}, Nothing, Nothing, Nothing}
const GLUVWMesh = UVWMesh{Float32, GLTriangle, Float32}

const NormalMesh{VT, FT, NT} = HMesh{Point{3, VT}, FT, Normal{3, NT}, Nothing, Nothing, Nothing, Nothing}
const GLNormalMesh = NormalMesh{Float32, GLTriangle, Float32}

const UVMesh2D{VT, FT, UVT} = HMesh{Point{2, VT}, FT, Nothing, UV{UVT}, Nothing, Nothing, Nothing}
const GLUVMesh2D = UVMesh2D{Float32, GLTriangle, Float32}

const NormalColorMesh{VT, FT, NT, CT} = HMesh{Point{3, VT}, FT, Normal{3, NT}, Nothing, CT, Nothing, Nothing}
const GLNormalColorMesh = NormalColorMesh{Float32, GLTriangle, Float32, RGBA{Float32}}

const NormalVertexcolorMesh{VT, FT, NT, CT} = HMesh{Point{3, VT}, FT, Normal{3, NT}, Nothing, Vector{CT}, Nothing, Nothing}
const GLNormalVertexcolorMesh = NormalVertexcolorMesh{Float32, GLTriangle, Float32, RGBA{Float32}}

const NormalAttributeMesh{VT, FT, NT, AT, A_ID_T} = HMesh{Point{3, VT}, FT, Normal{3, NT}, Nothing, Nothing, AT, A_ID_T}
const GLNormalAttributeMesh = NormalAttributeMesh{Float32, GLTriangle, Float32, Vector{RGBA{N0f8}}, Float32}

const NormalUVWMesh{VT, FT, NT, UVT} = HMesh{Point{3, VT}, FT, Normal{3, NT}, UVW{UVT}, Nothing, Nothing, Nothing}
const GLNormalUVWMesh = NormalUVWMesh{Float32, GLTriangle, Float32, Float32}

const NormalUVMesh{VT, FT, NT, UVT} = HMesh{Point{3, VT}, FT, Normal{3, NT}, UV{UVT}, Nothing, Nothing, Nothing}
const GLNormalUVMesh = NormalUVMesh{Float32, GLTriangle, Float32, Float32}

for T in (
        :SimpleMesh,
        :PlainMesh,
        :GLPlainMesh,
        :Mesh2D,
        :GLMesh2D,
        :UVMesh,
        :GLUVMesh,
        :UVWMesh,
        :GLUVWMesh,
        :NormalMesh,
        :GLNormalMesh,
        :UVMesh2D,
        :GLUVMesh2D,
        :NormalColorMesh,
        :GLNormalColorMesh,
        :NormalVertexcolorMesh,
        :GLNormalVertexcolorMesh,
        :NormalAttributeMesh,
        :GLNormalAttributeMesh,
        :NormalUVWMesh,
        :GLNormalUVWMesh,
        :NormalUVMesh,
        :GLNormalUVMesh,
    )
    @eval Base.show(io::IO, ::Type{<: $T}) = print(io, $(string(T)))
end
