
"""
A point type that holds additional metadata
"""
struct MetaPoint{N, T, Names, Types} <: AbstractPoint{N, T}
    point::Point{N, T}
    meta::NamedTuple{Names, Types}
end

function Base.getproperty(x::MetaPoint{Dim, T, Names, Types}, field::Symbol) where {Dim, T, Names, Types}
    field === :polygon && return getfield(x, :polygon)
    Base.sym_in(field, Names) && return getfield(getfield(x, :meta), field)
    error("Field $field not part of Element")
end

struct MetaPolygon{
        Dimension, T <: Real,
        P <: AbstractPolygon{Dimension, T},
        Names, Types
    } <: AbstractPolygon{Dimension, T}

    polygon::P
    meta::NamedTuple{Names, Types}
end

MetaPolygon(polygon::AbstractPolygon; meta...) = MetaPolygon(polygon, values(meta))
function MetaPolygon(polygon::AbstractVector; meta...)
    MetaPolygon(Polygon(polygon); meta...)
end
function MetaPolygon(exterior::L, interior::AbstractVector{L}; meta...) where L
    MetaPolygon(Polygon(exterior, interior); meta...)
end

function Tables.schema(::AbstractVector{MetaPolygon{Dim, T, P, Names, Types}}) where {Dim, T, P, Names, Types}
    Tables.Schema((:polygon, Names...), (P, Types...))
end
function Base.getproperty(x::MetaPolygon{Dim, T, P, Names, Types}, field::Symbol) where {Dim, T, P, Names, Types}
    field === :polygon && return getfield(x, :polygon)
    Base.sym_in(field, Names) && return getfield(getfield(x, :meta), field)
    error("Field $field not part of Element")
end

Tables.rows(x::MultiPolygon) = x.polygons
Tables.istable(::Type{<:MultiPolygon}) = true
Tables.istable(::Type{<:MultiPolygon{<:Tuple}}) = false
Tables.rowaccess(::Type{<:MultiPolygon}) = true
# Tables.columnaccess(::Type{<:MultiPolygon}) = true
