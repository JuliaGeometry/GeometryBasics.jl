"""
OffsetInteger type mainly for indexing.
* `O` - The offset relative to Julia arrays. This helps reduce copying when
communicating with 0-indexed systems such ad OpenGL.
"""
struct OffsetInteger{O, T <: Integer} <: Integer
    i::T

    OffsetInteger{O, T}(x::Integer) where {O, T <: Integer} = new{O, T}(T(O >= 0 ? x + O : x - (-O)))
end
raw(x::OffsetInteger) = x.i
raw(x::Integer) = x

import Base: +, -, abs, *, /, div, ==, <=, >=, show, to_index

function show(io::IO, oi::OffsetInteger{O, T}) where {O, T}
    print(io, "|$(raw(oi)) (indexes as $(O >= 0 ? raw(oi) - O : raw(oi) + -O))|")
end

Base.eltype(::Type{OffsetInteger{O, T}}) where {O, T} = T
Base.eltype(oi::OffsetInteger) = eltype(typeof(oi))

# constructors and conversion
OffsetInteger{O1, T1}(x::OffsetInteger{O2, T2}) where {O1, O2, T1 <: Integer, T2 <: Integer} = OffsetInteger{O1, T1}(T2(x))

OffsetInteger{O}(x::Integer) where {O} = OffsetInteger{O, eltype(x)}(x)
OffsetInteger{O}(x::OffsetInteger) where {O} = OffsetInteger{O, eltype(x)}(x)
(::Type{IT})(x::OffsetInteger{O, T}) where {IT <: Integer, O, T <: Integer} = IT(raw(x) + -O)

Base.@pure pure_max(x1, x2) = x1 > x2 ? x1 : x2
Base.promote_rule(::Type{T1}, ::Type{OffsetInteger{O, T2}}) where {T1 <: Integer, O, T2} = T1
Base.promote_rule(::Type{OffsetInteger{O1, T1}}, ::Type{OffsetInteger{O2, T2}}) where {O1, O2, T1, T2} = OffsetInteger{pure_max(O1, O2), promote_type(T1, T2)}

to_index(I::OffsetInteger) = raw(OneIndex(I))
to_index(I::OffsetInteger{0}) = raw(I)

# basic operators
for op in (:(-), :abs)
    @eval $(op)(x::T) where {T <: OffsetInteger} = T($(op)(x.i))
end
for op in (:(+), :(-), :(*), :(/), :div)
    @eval begin
        @inline function $(op)(x::OffsetInteger{O}, y::OffsetInteger{O}) where O
            OffsetInteger{O}($op(x.i, y.i))
        end
    end
end
for op in (:(==), :(>=), :(<=))
    @eval begin
        @inline function $(op)(x::OffsetInteger{O}, y::OffsetInteger{O}) where O
            $op(x.i, y.i)
        end
    end
end
