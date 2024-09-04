import LinearAlgebra: inv, det
import Random

struct Mat{Row, Column, T, L} <: AbstractMatrix{T}
    values::NTuple{L, T}
    function Mat{R, C, T}(values::NTuple{L, T}) where {R, C, T, L}
        @assert L == R * C "$R * $C needs to be $L"
        return new{R, C, T, L}(values)
    end
        function Mat{R, C, T, L}(values::NTuple{L, T}) where {R, C, T, L}
        @assert L == R * C "$R * $C needs to be $L"
        return new{R, C, T, L}(values)
    end
end

Base.size(::Mat{R, C}) where {R, C} = (R, C)
Base.size(::Type{<: Mat{R, C}}) where {R, C} = (R, C)
Base.ndims(::Type{<: Mat}) = 2
Base.getindex(mat::Mat{R, C}, i) where {R, C} = mat.values[i]

# TODO: maybe ranges as well?
function Base.getindex(mat::Mat{R, C}, i::StaticVector{N}, j::Integer) where {R, C, N}
    @boundscheck begin
        all(x -> 1 ≤ x ≤ R, i) && (1 ≤ j ≤ C)
    end
    @inbounds data = ntuple(n -> mat.values[i[n] + R * (j-1)], N)
    return Mat{length(i), 1}(data)
end
function Base.getindex(mat::Mat{R, C}, i::Integer, j::StaticVector{N}) where {R, C, N}
    @boundscheck begin
        (1 ≤ i ≤ R) && all(x -> 1 ≤ x ≤ C, j)
    end
    @inbounds data = ntuple(n -> mat.values[i + R * (j[n] - 1)], N)
    return Mat{1, length(j)}(data)
end
function Base.getindex(mat::Mat{R, C}, i::StaticVector{N}, j::StaticVector{M}) where {R, C, N, M}
    @boundscheck begin
        all(x -> 1 ≤ x ≤ R, i) && all(x -> 1 ≤ x ≤ C, j)
    end
    data = ntuple(N * M) do nm
        m, n = fldmod1(nm, R)
        @inbounds return mat.values[i[n] + R * (j[m] - 1)]
    end
    return Mat{1, length(j)}(data)
end

Base.IndexStyle(::Mat)= Base.IndexLinear()

function Mat{C, R, T}(::LinearAlgebra.UniformScaling) where {C, R, T}
    idx = CartesianIndices((R, C))
    data = ntuple(C * R) do i
        ci, ri = Tuple(idx[i])
        return ci === ri ? T(1) : T(0)
    end
    return Mat{R, C, T}(data)
end

function Mat{C, R, T, L}(::LinearAlgebra.UniformScaling) where {C, R, T, L}
    idx = CartesianIndices((R, C))
    data = ntuple(C * R) do i
        ci, ri = Tuple(idx[i])
        return ci === ri ? T(1) : T(0)
    end
    return Mat{R, C, T, L}(data)
end

Mat{C, R, T}(args...) where {C, R, T} = Mat{C, R, T}(args)
Mat{C, R, T, L}(args...) where {C, R, T, L} = Mat{C, R, T}(args)
Mat{C, R, T, L}(mat::Mat{C, R}) where {C, R, T, L} = Mat{C, R, T}(mat.values)
Mat{C}(args...) where {C} = Mat{C, C}(args)
Mat{C}(arg) where {C} = Mat{C, C}(arg)
Mat{C, R}(x::Tuple) where {C, R} = Mat{C, R}(promote(x...))
Mat{C, R, T}(x::Tuple) where {C, R, T} = Mat{C, R, T}(convert(NTuple{length(x), T}, x))
Mat{C, R}(x::NTuple{L, T}) where {C, R, L, T} = Mat{C, R, T}(x)
Mat{C, R, T1}(x::NTuple{L, T2}) where {C, R, L, T1, T2} = Mat{C, R, T1}(convert(NTuple{L, T1}, x))

# catch Mat2(...) etc which lowers to Mat{2, 2, T, 4} where T
(::Type{MT})(args...) where {C, R, MT <: Mat{C, R}} = Mat{C, R}(args)
(::Type{MT})(args::Tuple) where {C, R, MT <: Mat{C, R}} = Mat{C, R}(args)

Mat{C, R}(x::AbstractMatrix{T}) where {C, R, T} = Mat{C, R, T}(x)
Mat{C, R, T}(x::AbstractMatrix) where {C, R, T} = Mat{C, R, T}(ntuple(i-> convert(T, x[i]), C*R))

Base.convert(::Type{Mat{C, R, T, L}}, from::Mat{C, R}) where {C, R, T, L} = Mat{C, R, T}(from.values)

# Matrix products
# General shape mismatched versions are errors
(*)(a::Mat{M, N, T1}, b::Mat{O, K, T2}) where {T1, T2, M, N, O, K} = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))

Base.:(-)(a::Mat{R, C}, b::Mat{R, C}) where {R, C} = Mat{R, C}(a.values .- b.values)
Base.:(+)(a::Mat{R, C}, b::Mat{R, C}) where {R, C} = Mat{R, C}(a.values .+ b.values)

# matrix * matrix
@generated function *(a::Mat{M, N, T1}, b::Mat{N, P, T2}) where {T1, T2, M, N, P}
    elements = Expr(:tuple)
    for j in 1:P
        for i in 1:M
            plus_expr = Expr(:call, :+, (:(a[$i,$k]*b[$k,$j]) for k = 1:N)...)
            push!(elements.args, plus_expr)
        end
    end
    :(Mat{$M, $P}($elements))
end

@inline det(A::Mat{1, 1, T}) where {T} = @inbounds return ( A[1] )
@inline det(A::Mat{2, 2, T}) where {T} = @inbounds return ( A[1,1]*A[2,2] - A[1,2]*A[2,1])
@inline det(A::Mat{3, 3, T}) where {T} = @inbounds return (
    A[1,1]*(A[2,2]*A[3,3]-A[2,3]*A[3,2]) -
    A[1,2]*(A[2,1]*A[3,3]-A[2,3]*A[3,1]) +
    A[1,3]*(A[2,1]*A[3,2]-A[2,2]*A[3,1])
)

det(A::Mat{4, 4, T}) where {T} = @inbounds return (
    A[13] * A[10]  * A[7]  * A[4]  - A[9] * A[14] * A[7]  * A[4]   -
    A[13] * A[6]   * A[11] * A[4]  + A[5] * A[14] * A[11] * A[4]   +
    A[9]  * A[6]   * A[15] * A[4]  - A[5] * A[10] * A[15] * A[4]   -
    A[13] * A[10]  * A[3]  * A[8]  + A[9] * A[14] * A[3]  * A[8]   +
    A[13] * A[2]   * A[11] * A[8]  - A[1] * A[14] * A[11] * A[8]   -
    A[9]  * A[2]   * A[15] * A[8]  + A[1] * A[10] * A[15] * A[8]   +
    A[13] * A[6]   * A[3]  * A[12] - A[5] * A[14] * A[3]  * A[12]  -
    A[13] * A[2]   * A[7]  * A[12] + A[1] * A[14] * A[7]  * A[12]  +
    A[5]  * A[2]   * A[15] * A[12] - A[1] * A[6]  * A[15] * A[12]  -
    A[9]  * A[6]   * A[3]  * A[16] + A[5] * A[10] * A[3]  * A[16]  +
    A[9]  * A[2]   * A[7]  * A[16] - A[1] * A[10] * A[7]  * A[16]  -
    A[5]  * A[2]   * A[11] * A[16] + A[1] * A[6]  * A[11] * A[16]
)

det(A::Mat) = det(Matrix(A))

inv(A::Mat{1, 1, T, 1}) where T = @inbounds return Mat{1, 1, T, 1}(inv(A[1]))
function inv(A::Mat{2, 2, T, L}) where {T, L}
    determinant = det(A)
    @inbounds return Mat{2, 2, T}(
        (A[2,2] /determinant, -A[2,1]/determinant),
        (-A[1,2]/determinant, A[1,1] /determinant)
    )
end
function inv(A::Mat{3, 3, T, L}) where {T, L}
    determinant = det(A)
    @inbounds return  Mat{3, 3, T}(
        (A[2,2]*A[3,3]-A[2,3]*A[3,2]) /determinant,
        -(A[2,1]*A[3,3]-A[2,3]*A[3,1])/determinant,
        (A[2,1]*A[3,2]-A[2,2]*A[3,1]) /determinant,

        -(A[1,2]*A[3,3]-A[1,3]*A[3,2])/determinant,
        (A[1,1]*A[3,3]-A[1,3]*A[3,1]) /determinant,
        -(A[1,1]*A[3,2]-A[1,2]*A[3,1])/determinant,

        (A[1,2]*A[2,3]-A[1,3]*A[2,2]) /determinant,
        -(A[1,1]*A[2,3]-A[1,3]*A[2,1])/determinant,
        (A[1,1]*A[2,2]-A[1,2]*A[2,1]) /determinant
    )
end

function inv(A::Mat{4, 4, T, L}) where {T, L}
    determinant = det(A)
    @inbounds return Mat{4, 4, T}(
        (A[2,3]*A[3,4]*A[4,2] - A[2,4]*A[3,3]*A[4,2] + A[2,4]*A[3,2]*A[4,3] - A[2,2]*A[3,4]*A[4,3] - A[2,3]*A[3,2]*A[4,4] + A[2,2]*A[3,3]*A[4,4]) / determinant,
        (A[2,4]*A[3,3]*A[4,1] - A[2,3]*A[3,4]*A[4,1] - A[2,4]*A[3,1]*A[4,3] + A[2,1]*A[3,4]*A[4,3] + A[2,3]*A[3,1]*A[4,4] - A[2,1]*A[3,3]*A[4,4]) / determinant,
        (A[2,2]*A[3,4]*A[4,1] - A[2,4]*A[3,2]*A[4,1] + A[2,4]*A[3,1]*A[4,2] - A[2,1]*A[3,4]*A[4,2] - A[2,2]*A[3,1]*A[4,4] + A[2,1]*A[3,2]*A[4,4]) / determinant,
        (A[2,3]*A[3,2]*A[4,1] - A[2,2]*A[3,3]*A[4,1] - A[2,3]*A[3,1]*A[4,2] + A[2,1]*A[3,3]*A[4,2] + A[2,2]*A[3,1]*A[4,3] - A[2,1]*A[3,2]*A[4,3]) / determinant,

        (A[1,4]*A[3,3]*A[4,2] - A[1,3]*A[3,4]*A[4,2] - A[1,4]*A[3,2]*A[4,3] + A[1,2]*A[3,4]*A[4,3] + A[1,3]*A[3,2]*A[4,4] - A[1,2]*A[3,3]*A[4,4]) / determinant,
        (A[1,3]*A[3,4]*A[4,1] - A[1,4]*A[3,3]*A[4,1] + A[1,4]*A[3,1]*A[4,3] - A[1,1]*A[3,4]*A[4,3] - A[1,3]*A[3,1]*A[4,4] + A[1,1]*A[3,3]*A[4,4]) / determinant,
        (A[1,4]*A[3,2]*A[4,1] - A[1,2]*A[3,4]*A[4,1] - A[1,4]*A[3,1]*A[4,2] + A[1,1]*A[3,4]*A[4,2] + A[1,2]*A[3,1]*A[4,4] - A[1,1]*A[3,2]*A[4,4]) / determinant,
        (A[1,2]*A[3,3]*A[4,1] - A[1,3]*A[3,2]*A[4,1] + A[1,3]*A[3,1]*A[4,2] - A[1,1]*A[3,3]*A[4,2] - A[1,2]*A[3,1]*A[4,3] + A[1,1]*A[3,2]*A[4,3]) / determinant,

        (A[1,3]*A[2,4]*A[4,2] - A[1,4]*A[2,3]*A[4,2] + A[1,4]*A[2,2]*A[4,3] - A[1,2]*A[2,4]*A[4,3] - A[1,3]*A[2,2]*A[4,4] + A[1,2]*A[2,3]*A[4,4]) / determinant,
        (A[1,4]*A[2,3]*A[4,1] - A[1,3]*A[2,4]*A[4,1] - A[1,4]*A[2,1]*A[4,3] + A[1,1]*A[2,4]*A[4,3] + A[1,3]*A[2,1]*A[4,4] - A[1,1]*A[2,3]*A[4,4]) / determinant,
        (A[1,2]*A[2,4]*A[4,1] - A[1,4]*A[2,2]*A[4,1] + A[1,4]*A[2,1]*A[4,2] - A[1,1]*A[2,4]*A[4,2] - A[1,2]*A[2,1]*A[4,4] + A[1,1]*A[2,2]*A[4,4]) / determinant,
        (A[1,3]*A[2,2]*A[4,1] - A[1,2]*A[2,3]*A[4,1] - A[1,3]*A[2,1]*A[4,2] + A[1,1]*A[2,3]*A[4,2] + A[1,2]*A[2,1]*A[4,3] - A[1,1]*A[2,2]*A[4,3]) / determinant,

        (A[1,4]*A[2,3]*A[3,2] - A[1,3]*A[2,4]*A[3,2] - A[1,4]*A[2,2]*A[3,3] + A[1,2]*A[2,4]*A[3,3] + A[1,3]*A[2,2]*A[3,4] - A[1,2]*A[2,3]*A[3,4]) / determinant,
        (A[1,3]*A[2,4]*A[3,1] - A[1,4]*A[2,3]*A[3,1] + A[1,4]*A[2,1]*A[3,3] - A[1,1]*A[2,4]*A[3,3] - A[1,3]*A[2,1]*A[3,4] + A[1,1]*A[2,3]*A[3,4]) / determinant,
        (A[1,4]*A[2,2]*A[3,1] - A[1,2]*A[2,4]*A[3,1] - A[1,4]*A[2,1]*A[3,2] + A[1,1]*A[2,4]*A[3,2] + A[1,2]*A[2,1]*A[3,4] - A[1,1]*A[2,2]*A[3,4]) / determinant,
        (A[1,2]*A[2,3]*A[3,1] - A[1,3]*A[2,2]*A[3,1] + A[1,3]*A[2,1]*A[3,2] - A[1,1]*A[2,3]*A[3,2] - A[1,2]*A[2,1]*A[3,3] + A[1,1]*A[2,2]*A[3,3]) / determinant
    )
end

inv(A::Mat) = typeof(A)(inv(Matrix(A)))

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<: Mat{R, C, T}}) where {R,C,T}
    return Mat{R, C, T}(ntuple(i-> rand(rng, T), R*C))
end

function Random.randn(rng::Random.AbstractRNG, ::Type{<: Mat{R, C, T}}) where {R,C,T}
    return Mat{R, C, T}(ntuple(i-> randn(rng, T), R*C))
end

function Base.transpose(a::Mat{R, C, T}) where {R, C, T}
    idx = CartesianIndices((R, C))
    data = ntuple(R * C) do i
        ci, ri = Tuple(idx[i])
        return a[ri, ci]
    end
    Mat{R, C, T}(data)
end

(*)(a::Mat{M, N, T1}, b::StaticVector{O, T2}) where {T1, T2, M, N, O} = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))

# vector * (row vector)
@generated function *(a::StaticVector{N, T1}, b::Mat{1, M, T2}) where {T1, T2, N, M}
    elements = Expr(:tuple, [Expr(:tuple, [:(a[$i] * b[$j]) for i in 1:N]...) for j in 1:M]...)
    return :(similar_type(a)($elements))
end

# matrix * vector
@generated function *(a::Mat{M, N, T1}, b::StaticVector{N, T2}) where {T1, T2, M, N}
    elements = Expr(:tuple, [Expr(:call, :+, [:(a[$i,$k]*b[$k]) for k = 1:N]...) for i in 1:M]...)
    return :(similar_type(b)($elements))
end

# TODO: delete since we have a more specialized version up top?
function Base.getindex(mat::Mat{R, C, T}, r::Vec{NR}, c::Vec{NC}) where {R, C, NR, NC, T}
    idx = CartesianIndices((NR, NC))
    data = ntuple(NR * NC) do i
        ri, ci = Tuple(idx[i])
        return mat[r[ri], c[ci]]
    end
    return Mat{NR, NC, T}(data)
end

# TODO: Fix Vec(mat) becoming Vec((mat,)) (i.e. restrict eltype to Number?)
(VT::Type{<: StaticVector{N}})(mat::Mat{N, 1}) where {N} = VT(mat.values)

function Base.isapprox(
        a::Mat{R1, C1, T1}, b::Mat{R2, C2, T2};
        atol::Real = 0,
        rtol::Real = atol > 0 ? 0 : sqrt(max(eps(T1), eps(T2)))
    ) where {R1, R2, C1, C2, T1, T2}
    return (R1 == R2) && (C1 == C2) &&
        norm(a - b) <= max(atol, rtol * max(norm(a), norm(b)))
end