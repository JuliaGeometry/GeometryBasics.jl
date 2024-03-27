
function unit(::Type{T}, i::Integer) where {T <: StaticVector}
    tup = ntuple(Val(length(T))) do j
        return ifelse(i == j, 1, 0)
    end
    return T(tup)
end

macro fixed_vector(name_parent)
    @assert name_parent.head == :(=)
    name, parent = name_parent.args
    expr = quote
        struct $(name){S,T} <: $(parent){S,T}
            data::NTuple{S,T}

            function $(name){S,T}(x::NTuple{S,T}) where {S,T}
                return new{S,T}(x)
            end

            function $(name){S,T}(x::NTuple{S,Any}) where {S,T}
                return new{S,T}(StaticArrays.convert_ntuple(T, x))
            end
        end

        size_or(::Type{$(name)}, or) = or
        eltype_or(::Type{$(name)}, or) = or
        eltype_or(::Type{$(name){S,T} where S}, or) where {T} = T
        eltype_or(::Type{$(name){S,T} where T}, or) where {S} = or
        eltype_or(::Type{$(name){S,T}}, or) where {S,T} = T

        size_or(::Type{$(name){S,T} where S}, or) where {T} = or
        size_or(::Type{$(name){S,T} where T}, or) where {S} = Size{(S,)}()
        size_or(::Type{$(name){S,T}}, or) where {S,T} = (S,)

        # Array constructor
        function $(name){S}(x::AbstractVector{T}) where {S,T}
            @assert S <= length(x)
            return $(name){S,T}(ntuple(i -> x[i], Val(S)))
        end

        function $(name){S,T1}(x::AbstractVector{T2}) where {S,T1,T2}
            @assert S <= length(x)
            return $(name){S,T1}(ntuple(i -> convert(T1, x[i]), Val(S)))
        end

        function $(name){S,T}(x) where {S,T}
            return $(name){S,T}(ntuple(i -> convert(T, x), Val(S)))
        end

        $(name){S}(x::T) where {S,T} = $(name){S,T}(ntuple(i -> x, Val(S)))
        $(name){1,T}(x::T) where {T} = $(name){1,T}((x,))
        $(name)(x::NTuple{S}) where {S} = $(name){S}(x)
        function $(name)(x::T) where {S,T <: Tuple{Vararg{Any,S}}}
            return $(name){S,StaticArrays.promote_tuple_eltype(T)}(x)
        end

        function $(name){S}(x::T) where {S,T <: Tuple}
            return $(name){S,StaticArrays.promote_tuple_eltype(T)}(x)
        end
        $(name){S,T}(x::StaticVector) where {S,T} = $(name){S,T}(Tuple(x))

        @generated function (::Type{$(name){S,T}})(x::$(name)) where {S,T}
            idx = [:(x[$i]) for i in 1:S]
            return quote
                $($(name)){S,T}($(idx...))
            end
        end

        @generated function Base.convert(::Type{$(name){S,T}}, x::$(name)) where {S,T}
            idx = [:(x[$i]) for i in 1:S]
            return quote
                $($(name)){S,T}($(idx...))
            end
        end

        @generated function (::Type{SV})(x::StaticVector) where {SV <: $(name)}
            len = size_or(SV, size(x))[1]
            return if length(x) == len
                :(SV(Tuple(x)))
            elseif length(x) > len
                elems = [:(x[$i]) for i in 1:len]
                :(SV($(Expr(:tuple, elems...))))
            else
                error("Static Vector too short: $x, target type: $SV")
            end
        end

        Base.@pure StaticArrays.Size(::Type{$(name){S,Any}}) where {S} = Size(S)
        Base.@pure StaticArrays.Size(::Type{$(name){S,T}}) where {S,T} = Size(S)

        Base.@propagate_inbounds function Base.getindex(v::$(name){S,T}, i::Int) where {S,T}
            return v.data[i]
        end

        Base.Tuple(v::$(name)) = v.data
        function Base.convert(::Type{$(name){S,T}}, x::NTuple{S,T}) where {S,T}
            return $(name){S,T}(x)
        end
        function Base.convert(::Type{$(name){S,T}}, x::Tuple) where {S,T}
            return $(name){S,T}(convert(NTuple{S,T}, x))
        end

        @generated function StaticArrays.similar_type(::Type{SV}, ::Type{T},
                                                      s::Size{S}) where {SV <: $(name),T,S}
            return if length(S) === 1
                $(name){S[1],T}
            else
                StaticArrays.default_similar_type(T, s(), Val{length(S)})
            end
        end

        Base.:(*)(a::$name, b::$name) = a .* b
        function Base.broadcasted(f, a::AbstractArray{T}, b::$name) where {T <: $name}
            return Base.broadcasted(f, a, (b,))
        end
    end
    return esc(expr)
end

abstract type AbstractPoint{Dim,T} <: StaticVector{Dim,T} end

@fixed_vector Point = AbstractPoint
@fixed_vector Vec = StaticVector

Base.isnan(p::Union{AbstractPoint,Vec}) = any(x -> isnan(x), p)

const Mat = SMatrix
const VecTypes{N,T} = Union{StaticVector{N,T},NTuple{N,T}}
const Vecf{N} = Vec{N,Float32}
const Pointf{N} = Point{N,Float32}

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

include("generated-aliases.jl")

export Mat, Vec, Point, unit
export Vecf, Pointf
