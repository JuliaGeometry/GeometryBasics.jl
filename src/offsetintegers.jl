
"""
OffsetInteger type mainly for indexing.
* `O` - The offset relative to Julia arrays. This helps reduce copying when
communicating with 0-indexed systems such ad OpenGL.
"""
struct OffsetInteger{O, T <: Integer} <: Integer
    i::T
    OffsetInteger{O, T}(x::Integer) where {O, T <: Integer} = new{O, T}(T(x + O))
end

const ZeroIndex{T <: Integer} = OffsetInteger{-1, T}
const OneIndex{T <: Integer} = OffsetInteger{0, T}
const GLIndex = ZeroIndex{Cuint}

raw(x::OffsetInteger) = x.i
raw(x::Integer) = x
value(x::OffsetInteger{O, T}) where {O, T} = raw(x) - O
value(x::Integer) = x

function show(io::IO, oi::OffsetInteger{O, T}) where {O, T}
    print(io, "|$(raw(oi)) (indexes as $(value(oi))|")
end

Base.eltype(::Type{OffsetInteger{O, T}}) where {O, T} = T
Base.eltype(oi::OffsetInteger) = eltype(typeof(oi))

# constructors and conversion
OffsetInteger{O1, T1}(x::OffsetInteger{O2, T2}) where {O1, O2, T1 <: Integer, T2 <: Integer} = OffsetInteger{O1, T1}(T2(x))

OffsetInteger{O}(x::Integer) where {O} = OffsetInteger{O, eltype(x)}(x)
OffsetInteger{O}(x::OffsetInteger) where {O} = OffsetInteger{O, eltype(x)}(x)
(::Type{IT})(x::OffsetInteger{O, T}) where {IT <: Integer, O, T <: Integer} = IT(value(x))

Base.@pure pure_max(x1, x2) = x1 > x2 ? x1 : x2
Base.promote_rule(::Type{T1}, ::Type{OffsetInteger{O, T2}}) where {T1 <: Integer, O, T2} = T1
Base.promote_rule(::Type{OffsetInteger{O1, T1}}, ::Type{OffsetInteger{O2, T2}}) where {O1, O2, T1, T2} = OffsetInteger{pure_max(O1, O2), promote_type(T1, T2)}

#Need to convert to Int here because of: https://github.com/JuliaLang/julia/issues/35038
Base.to_index(I::OffsetInteger) = convert(Int, raw(OneIndex(I)))
Base.to_index(I::OffsetInteger{0}) = convert(Int, raw(I))

# basic operators
for op in (:(-), :abs)
    @eval Base.$op(x::T) where {T <: OffsetInteger} = T($(op)(value(x)))
end

for op in (:(+), :(-), :(*), :(/), :div)
    @eval begin
        @inline function Base.$op(x::OffsetInteger{O}, y::OffsetInteger{O}) where O
            OffsetInteger{O}($op(value(x), value(y)))
        end
    end
end

for op in (:(==), :(>=), :(<=), :(<) , :(>), :sub_with_overflow)
    @eval begin
        @inline function Base.$op(x::OffsetInteger{O}, y::OffsetInteger{O}) where O
            $op(x.i, y.i)
        end
    end
end
