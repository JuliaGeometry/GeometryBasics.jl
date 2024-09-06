
using LinearAlgebra
import Random
import Base: setindex

abstract type StaticVector{N, T} <: AbstractVector{T} end
function similar_type end

struct StaticArrayStyle{T, AsConst} <: Broadcast.BroadcastStyle end
StaticArrayStyle{T}() where T = StaticArrayStyle{T, false}()

macro fixed_vector(name_parent)
    @assert name_parent.head == :(=)
    VecT, SuperT = name_parent.args

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

        $(VecT){1, T}(x::Number) where {T} = $(VecT){1, T}((x,))
        $(VecT){1, T}(x::Tuple{Any}) where T = $(VecT){1, T}((T(x[1]),))
        $(VecT)(x::Tuple) = $(VecT)(promote(x...))
        $(VecT){N}(x::Tuple) where {N} = $(VecT){N}(promote(x...))
        $(VecT){N, T}(x::Tuple) where {N,T} = $(VecT){N,T}(convert(NTuple{N,T}, x))

        $(VecT)(x::NTuple{N, T}) where {N, T} = $(VecT){N,T}(x)
        $(VecT){N}(x::NTuple{N, T}) where {N, T} = $(VecT){N,T}(x)

        $(VecT)(x::Vararg{Any,N}) where {N} = $(VecT){N}(x)
        $(VecT)(x::Vararg{T,N}) where {T,N} = $(VecT){N,T}(x)

        $(VecT){N}(x::Vararg{Any,N}) where {N} = $(VecT){N}(x)
        $(VecT){N}(x::Vararg{T,N}) where {T,N} = $(VecT){N,T}(x)

        $(VecT){N, T}(x::Vararg{Any,N}) where {T,N} = $(VecT){N,T}(x)
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

        @inline similar_type(::$(VecT){N, T}, n::Integer) where {N, T} = $(VecT){n}
        @inline similar_type(::$(VecT){N}, ::Type{T}) where {N, T} = $(VecT){N, T}
        @inline similar_type(::$(VecT), n::Integer, ::Type{T}) where {T} = $(VecT){n, T}
        @inline similar_type(::$(VecT)) = $(VecT)

        Base.BroadcastStyle(::Type{<: $(VecT)}) = StaticArrayStyle{$(VecT)}()
        Base.values(v::$(VecT)) = v.data

        function LinearAlgebra.cross(a::$(VecT){3}, b::$(VecT){3})
            @inbounds elements = (a[2]*b[3]-a[3]*b[2],
                                    a[3]*b[1]-a[1]*b[3],
                                    a[1]*b[2]-a[2]*b[1])
            return $(VecT)(elements)
        end
    end
    return esc(expr)
end

# Broadcasting
# style rules
Base.BroadcastStyle(::StaticArrayStyle{T, B1}, ::StaticArrayStyle{T, B2}) where {B1, B2, T} = StaticArrayStyle{T, B1 || B2}()
Base.BroadcastStyle(s::StaticArrayStyle, ::Broadcast.AbstractArrayStyle{0}) = s
Base.BroadcastStyle(::Broadcast.AbstractArrayStyle{0}, s::StaticArrayStyle) = s
Base.BroadcastStyle(s::StaticArrayStyle, ::Broadcast.Style{Tuple}) = s
Base.BroadcastStyle(::Broadcast.Style{Tuple}, s::StaticArrayStyle) = s
Base.BroadcastStyle(::StaticArrayStyle{T, B}, ::Broadcast.BroadcastStyle) where {B, T} = StaticArrayStyle{T, true}()
Base.BroadcastStyle(::Broadcast.BroadcastStyle, ::StaticArrayStyle{T, B}) where {B, T} = StaticArrayStyle{T, true}()
# to allow mixing types, define:
# Base.BroadcastStyle(::StaticArrayStyle{<: Type1, B1}, ::StaticArrayStyle{<: Type2, B2}) where {B1, B2} = 
#     StaticArrayStyle{preffered_type, B1 || B2}()
# Base.BroadcastStyle(::StaticArrayStyle{<: Type2, B1}, ::StaticArrayStyle{<: Type1, B2}) where {B1, B2} = 
#     StaticArrayStyle{preffered_type, B1 || B2}()

# If we don't inherit from AbstractVector we need this?
# Base.broadcastable(x::StaticVector) = x

# Required to avoid size missmatches between Array and StaticVector
function Broadcast.instantiate(bc::Broadcast.Broadcasted{<: StaticArrayStyle{<: Any, true}})
    # transform this to an Array broadcast with Ref'd StaticVectors and tuples
    args_converted = map(arg -> arg isa Broadcast.Broadcasted ? copy(Broadcast.instantiate(arg)) : arg, bc.args)
    maybe_const_args = map(args_converted) do arg
        style = Base.BroadcastStyle(typeof(arg))
        if style isa Broadcast.AbstractArrayStyle # value or Array
            return arg
        else # tuple, StaticVector
            return Ref(arg)
        end
    end
    return Broadcast.broadcasted(bc.f, maybe_const_args...)
end

# resolve element-wise operation
function Base.copy(bc::Broadcast.Broadcasted{StaticArrayStyle{T, false}}) where T
    # Broadcasted may end up in args from nested calls (e.g. foo(a, b .+ c); a .+ b .+ c)
    args = map(arg -> values(arg isa Broadcast.Broadcasted ? copy(arg) : arg), bc.args)
    return T(broadcast(bc.f, args...))
end

Base.map(f, a::StaticVector, args::AbstractArray...) = broadcast(f, a, args...)
Base.map(f, a::AbstractArray, b::StaticVector, args::AbstractArray...) = broadcast(f, a, b, args...)
Base.map(f, a::StaticVector, b::StaticVector, args::AbstractArray...) = broadcast(f, a, b, args...)

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{V}) where V <: StaticVector{N,T} where {N, T}
    V(ntuple(i-> rand(rng, T), N))
end
function Random.randn(rng::Random.AbstractRNG, ::Type{V}) where V <: StaticVector{N,T} where {N, T}
    V(ntuple(i-> randn(rng, T), N))
end

Base.@propagate_inbounds function Base.getindex(v::StaticVector{N,T}, i::Integer) where {N,T}
    return v.data[i]
end

Base.setindex(c::V, v, i::Integer) where {V <: StaticVector} = V(Base.setindex(c.data, v, i))

Base.@propagate_inbounds function Base.getindex(a::StaticVector{N1, T}, idx::StaticVector{N, <:Integer}) where {N,N1,T}
    return similar_type(idx, N, T)(map(i-> a[i], idx))
end

Base.:(-)(a::StaticVector) = (-).(a)

import Base: *, +, -, /

for op in [:*, :+, :-, :/]
    @eval begin
        ($op)(a::StaticVector, b::StaticVector) = Base.broadcast($(op), a, b)
        ($op)(a::Number, b::StaticVector) = Base.broadcast($(op), a, b)
        ($op)(a::StaticVector, b::Number) = Base.broadcast($(op), a, b)
    end
end

LinearAlgebra.cross(a::StaticVector{2}, b::StaticVector{2}) = a[1]*b[2]-a[2]*b[1]
LinearAlgebra.norm(a::StaticVector) = sqrt(dot(a,a))
LinearAlgebra.normalize(a::StaticVector) = a ./ norm(a)

Base.eltype(::StaticVector{N, T}) where {N, T} = T
Base.eltype(::Type{<: StaticVector{N, T}}) where {N, T} = T

Base.size(::StaticVector{N}) where {N} = (N,)
Base.size(::Type{<: StaticVector{N}}) where {N} = (N,)
Base.length(::StaticVector{N}) where {N} = N
Base.length(::Type{<: StaticVector{N}}) where {N} = N
Base.ndims(::Type{<: StaticVector}) = 1

Base.copy(v::StaticVector) = deepcopy(v)

Base.vcat(a::StaticVector, b::StaticVector) = (a..., b...)

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
Base.zero(::V) where {V <:StaticVector} = zero(V)

function Base.:(==)(a::StaticVector{N}, b::StaticVector{N}) where N
    for i in 1:N
        a[i] == b[i] || return false
    end
    return true
end

function Base.isapprox(
        a::StaticVector{N1, T1}, b::StaticVector{N2, T2};
        atol::Real = 0,
        rtol::Real = atol > 0 ? 0 : sqrt(max(eps(T1), eps(T2)))
    ) where {N1, N2, T1, T2}
    return (N1 == N2) && norm(a - b) <= max(atol, rtol * max(norm(a), norm(b)))
end

@generated function Base.transpose(b::StaticVector{N,T}) where {N,T}
    expr = [:(transpose(b[$i])) for i=1:N]
    return quote
        Mat{1,N,T}($(expr...))
    end
end
Base.reverse(x::P) where P <: StaticVector = P(reverse(x.data))

# Since we don't inherit from AbstractArray, some extra functions need to be overloaded
LinearAlgebra.promote_leaf_eltypes(x::StaticVector{N, T}) where {N,T} = T

@fixed_vector Point = StaticVector
@fixed_vector Vec = StaticVector


Base.lastindex(::StaticVector{N}) where N = N

# Allow mixing Point Vec in broadcast
Base.BroadcastStyle(::StaticArrayStyle{<: Point, B1}, ::StaticArrayStyle{<: Vec, B2}) where {B1, B2} = 
    StaticArrayStyle{Point, B1 || B2}()
Base.BroadcastStyle(::StaticArrayStyle{<: Vec, B1}, ::StaticArrayStyle{<: Point, B2}) where {B1, B2} = 
    StaticArrayStyle{Point, B1 || B2}()

Base.:(+)(a::Vec{N}, b::Point{N}) where {N} = Point{N}(a.data .+ b.data)

const VecTypes{N,T} = Union{StaticVector{N,T}, NTuple{N,T}}
const Vecf{N} = Vec{N, Float32}
const PointT{T} = Point{N,T} where N
const Pointf{N} = Point{N,Float32}
    
Base.isnan(p::Union{Point,Vec}) = any(isnan, p)
Base.isinf(p::Union{Point,Vec}) = any(isinf, p)
Base.isfinite(p::Union{Point,Vec}) = all(isfinite, p)

## Generate aliases
## As a text file instead of eval/macro, to not confuse code linter

#=
open(joinpath(@__DIR__, "generated-aliases.jl"), "w") do io
    for i in 1:4
        for T in [:Point, :Vec, :Mat]
            namei = "$T$i"
            res = T == :Mat ? "Mat{$i,$i,T,$(i * i)}" : "$T{$i,T}"
            println(io, "const $(namei){T} = $res")
            println(io, "export $namei")
            for (postfix, t) in ["d" => Float64, "f" => Float32, "i" => Int, "ui" => UInt]
                namep = "$T$i$postfix"
                println(io, "const $(namep) = $(namei){$t}")
                println(io, "export $namep")
                # mnamep = "$(mname)$postfix"
                # println(io, "const $mnamep = $mname{$t}")
                # println(io, "export $mnamep")
            end
        end
    end
end
=#

include("mat.jl")
include("generated-aliases.jl")

export Mat, Vec, Point, unit
export Vecf, Pointf
