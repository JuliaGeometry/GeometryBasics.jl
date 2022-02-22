
using LinearAlgebra
import Random
import Base: setindex

abstract type StaticVector{N, T} end

macro fixed_vector(VecT, SuperT)
    expr = quote
        struct $(VecT){N, T} <: $(SuperT){N, T}
            data::NTuple{N,T}

            function $(VecT){N,T}(x::NTuple{N,T}) where {N,T}
                return new{N,T}(x)
            end

            function $(VecT){N,T}(x::NTuple{N,Any}) where {N,T}
                return new{N,T}(convert(NTuple{N, T}, x))
            end
            $(VecT){1, T}(x::Tuple{T}) where {T} = new{1, T}(x)
        end

        # Array constructor
        $(VecT)(x::AbstractVector) = error("You need to supply size of vector")
        $(VecT){N}(x::AbstractVector{T}) where {N,T} = $(VecT){N,T}(x)
        function $(VecT){N,T1}(x::AbstractVector{T2}) where {N,T1,T2}
            @assert N <= length(x)
            return $(VecT){N,T1}(ntuple(i -> convert(T1, x[i]), N))
        end

        # StaticVector conversion
        $(VecT)(x::StaticVector{N, T}) where {N,T} = $(VecT){N, T}(x)
        $(VecT){N}(x::StaticVector{N2, T}) where {N,N2,T} = $(VecT){N,T}(x)
        function $(VecT){N1,T1}(x::StaticVector{N2, T2}) where {N1,T1,N2,T2}
            @assert N1 <= N2
            return $(VecT){N1,T1}(ntuple(i -> convert(T1, x[i]), N1))
        end
        function $(VecT){1,T1}(x::StaticVector{N2, T2}) where {T1,N2,T2}
            @assert 1 <= N2
            return $(VecT){1,T1}(ntuple(i -> convert(T1, x[i]), 1))
        end

        # repeat
        $(VecT){N}(x::T) where {N,T<:Number} = $(VecT){N, T}(x)
        function $(VecT){N,T}(x::Number) where {N,T}
            return $(VecT){N,T}(ntuple(i -> convert(T, x), N))
        end

        $(VecT){1, T}(x) where {T} = $(VecT){1, T}((x,))
        $(VecT){1, T}(x::Tuple{Any}) where T = $(VecT){1, T}((T(x[1]),))
        $(VecT)(x::Tuple) = $(VecT)(promote(x...))
        $(VecT){N}(x::Tuple) where {N} = $(VecT){N}(promote(x...))
        $(VecT){N, T}(x::Tuple) where {N,T} = $(VecT){N,T}(convert(NTuple{N,T}, x))

        $(VecT)(x::NTuple{N, T}) where {N, T} = $(VecT){N,T}(x)
        $(VecT){N}(x::NTuple{N, T}) where {N, T} = $(VecT){N,T}(x)

        $(VecT)(x::Vararg{<:Any,N}) where {N} = $(VecT){N}(x)
        $(VecT)(x::Vararg{T,N}) where {T,N} = $(VecT){N,T}(x)

        $(VecT){N}(x::Vararg{<:Any,N}) where {N} = $(VecT){N}(x)
        $(VecT){N}(x::Vararg{T,N}) where {T,N} = $(VecT){N,T}(x)

        $(VecT){N, T}(x::Vararg{<:Any,N}) where {T,N} = $(VecT){N,T}(x)
        $(VecT){N, T1}(x::Vararg{T2,N}) where {T1,T2,N} = $(VecT){N, T1}(x)

        Base.convert(::Type{$(VecT){N,T}}, x) where {N,T} = $(VecT){N,T}(x)
        Base.convert(::Type{$(VecT){N}}, x) where {N} = $(VecT){N}(x)
        Base.convert(::Type{$(VecT){N}}, x::$(VecT){N}) where {N} = x
        Base.convert(::Type{$(VecT){N,T}}, x::$(VecT){N,T}) where {N,T} = x

        function Base.convert(::Type{$(VecT){N,T}}, x::NTuple{N,T}) where {N,T}
            return $(VecT){N,T}(x)
        end
        function Base.convert(::Type{$(VecT){N,T}}, x::Tuple) where {N,T}
            return $(VecT){N,T}(convert(NTuple{N,T}, x))
        end
        Base.convert(::Type{$(VecT)}, x::Tuple) = $(VecT)(x)

        Base.@propagate_inbounds function Base.getindex(v::$(VecT){N,T}, i::Int) where {N,T}
            return v.data[i]
        end
        Base.setindex(c::$(VecT){N, T}, v, i::Integer) where {N,T} = $(VecT){N,T}(Base.setindex(c.data, v, i))
        Base.@propagate_inbounds function Base.getindex(a::AbstractArray{T}, idx::$(VecT){N, <:Integer}) where {N,T}
            return $(VecT){N,T}(map(i-> a[i], idx))
        end

        Base.@propagate_inbounds function Base.getindex(a::StaticVector{N1, T}, idx::$(VecT){N, <:Integer}) where {N,N1,T}
            return $(VecT){N,T}(map(i-> a[i], idx))
        end

        Base.Tuple(v::$(VecT)) = v.data


        function Base.broadcasted(f, a::AbstractArray{T}, b::$(VecT)) where {T <: $(VecT)}
            return Base.broadcasted(f, a, (b,))
        end

        function Base.broadcasted(f, a::$(VecT), b::$(VecT))
            return $(VecT)(map(f, a.data, b.data))
        end

        Base.broadcasted(f, a::$(VecT)) = $(VecT)(f.(a.data))
        Base.broadcasted(f, a::$(VecT), b) = $(VecT)(f.(a.data, b))
        Base.broadcasted(f, a, b::$(VecT)) = $(VecT)(f.(a, b.data))

        Base.map(f, b::$(VecT)) = $(VecT)(map(f, b.data))

        (*)(a::Mat{M, N, T1}, b::$(VecT){O, T2}) where {T1, T2, M, N, O} = throw(DimensionMismatch("$N != $O in $(typeof(a)) and $(typeof(b))"))

        # vector * (row vector)
        @generated function *(a::$VecT{N, T1}, b::Mat{1, M, T2}) where {T1, T2, N, M}
            elements = Expr(:tuple, [Expr(:tuple, [:(a[$i] * b[$j]) for i in 1:N]...) for j in 1:M]...)
            return :($($(VecT))($elements))
        end

        # matrix * vector
        @generated function *(a::Mat{M, N, T1}, b::$VecT{N, T2}) where {T1, T2, M, N}
            total_terms = M*N
            if total_terms <= 64
                # Full unrolling
                elements = Expr(:tuple, [Expr(:call, :+, [:(a[$i,$k]*b[$k]) for k = 1:N]...) for i in 1:M]...)
            else
                # Expand as a bunch of dot products
                elements = Expr(:tuple, [:(bilindot($($(VecT))(row(a,$i)),b)) for i in 1:M]...)
            end
            return :($($(VecT))($elements))
        end

        Base.:(*)(a::$VecT, b::$VecT) = a .* b
        Base.:(*)(a::Number, b::$VecT) = a .* b
        Base.:(*)(a::$VecT, b::Number) = a .* b

        Base.:(+)(a::$VecT, b::$VecT) = a .+ b
        Base.:(+)(a::Number, b::$VecT) = a .+ b
        Base.:(+)(a::$VecT, b::Number) = a .+ b

        Base.:(-)(a::$VecT) = (-).(a)
        Base.:(-)(a::$VecT, b::$VecT) = a .- b
        Base.:(-)(a::Number, b::$VecT) = a .- b
        Base.:(-)(a::$VecT, b::Number) = a .- b

        function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{$(VecT){N,T}}) where {N,T}
            $(VecT){N,T}(ntuple(i-> rand(rng, T), N))
        end
        function Random.randn(rng::Random.AbstractRNG, ::Type{$(VecT){N,T}}) where {N,T}
            $(VecT){N,T}(ntuple(i-> randn(rng, T), N))
        end
        function LinearAlgebra.cross(a::$(VecT){3}, b::$(VecT){3})
            @inbounds elements = (a[2]*b[3]-a[3]*b[2],
                                    a[3]*b[1]-a[1]*b[3],
                                    a[1]*b[2]-a[2]*b[1])
            return $(VecT)(elements)
        end
    end
    return esc(expr)
end

LinearAlgebra.cross(a::StaticVector{2}, b::StaticVector{2}) = a[1]*b[2]-a[2]*b[1]
LinearAlgebra.norm(a::StaticVector) = sqrt(dot(a,a))
LinearAlgebra.normalize(a::StaticVector) = a ./ norm(a)

Base.eltype(::StaticVector{N, T}) where {N, T} = T
Base.eltype(::Type{<: StaticVector{N, T}}) where {N, T} = T

Base.size(::StaticVector{N}) where {N} = (N,)
Base.length(::StaticVector{N}) where {N} = N
Base.length(::Type{<: StaticVector{N}}) where {N} = N
function Base.iterate(A::StaticVector, i=1)
    i - 1 < length(A) ? (A[i], i + 1) : nothing
end

function unit(::Type{T}, i::Integer) where {T <: StaticVector}
    tup = ntuple(Val(length(T))) do j
        return ifelse(i == j, 1, 0)
    end
    return T(tup)
end

Base.zero(::Type{V}) where {V <:StaticVector} = V(0)

function Base.:(==)(a::StaticVector{N}, b::StaticVector{N}) where N
    for i in 1:N
        a[i] == b[i] || return false
    end
    return true
end

function Base.isapprox(a::StaticVector{N}, b::StaticVector{N}; kw...) where N
    for i in 1:N
        isapprox(a[i], b[i]; kw...) || return false
    end
    return true
end

@generated function Base.transpose(b::StaticVector{N,T}) where {N,T}
    expr = [:(transpose(b[$i])) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end

# Since we don't inherit from AbstractArray, some extra functions need to be overloaded
LinearAlgebra.promote_leaf_eltypes(x::StaticVector{N, T}) where {N,T} = T

abstract type AbstractPoint{Dim,T} <: StaticVector{Dim,T} end

@fixed_vector Point AbstractPoint
@fixed_vector Vec StaticVector

Base.lastindex(::StaticVector{N}) where N = N

Base.broadcasted(f, a::Point, b::GeometryBasics.Vec) = Vec(f.(a.data, b.data))
Base.broadcasted(f, a::Vec, b::Point) = Vec(f.(a.data, b.data))

Base.:(+)(a::Vec{N}, b::Point{N}) where {N} = Point{N}(a.data .+ b.data)

const VecTypes{N,T} = Union{StaticVector{N,T}, NTuple{N,T}}
const Vecf{N} = Vec{N, Float32}
const Pointf{N} = Point{N,Float32}
Base.isnan(p::Union{AbstractPoint,Vec}) = any(x -> isnan(x), p)

for i in 1:4
    for T in [:Point, :Vec]
        name = Symbol("$T$i")
        namef = Symbol("$T$(i)f")
        @eval begin
            const $name = $T{$i}
            const $namef = $T{$i,Float32}
            export $name
            export $namef
        end
    end
    name = Symbol("Mat$i")
    namef = Symbol("Mat$(i)f")
    namef = Symbol("Mat$(i)f")
    @eval begin
        const $name{T} = $Mat{$i,$i,T, $(i*i)}
        const $namef = $name{Float32}
        export $name
        export $namef
    end
end

export Mat, Vec, Point, unit
export Vecf, Pointf
