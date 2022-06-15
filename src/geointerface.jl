# Implementation of trait based interface from https://github.com/JuliaGeo/GeoInterface.jl/

GeoInterface.isgeometry(::Type{<:AbstractGeometry}) = true
GeoInterface.isgeometry(::Type{<:AbstractFace}) = true
GeoInterface.isgeometry(::Type{<:AbstractPoint}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractGeometry}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractPoint}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:LineString}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractPolygon}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractFace}}) = true
GeoInterface.isgeometry(::Type{<:Mesh}) = true

GeoInterface.geomtrait(::Point) = GeoInterface.PointTrait()
GeoInterface.geomtrait(::Line) = GeoInterface.LineTrait()
GeoInterface.geomtrait(::LineString) = GeoInterface.LineStringTrait()
GeoInterface.geomtrait(::Polygon) = GeoInterface.PolygonTrait()
GeoInterface.geomtrait(::MultiPoint) = GeoInterface.MultiPointTrait()
GeoInterface.geomtrait(::MultiLineString) = GeoInterface.MultiLineStringTrait()
GeoInterface.geomtrait(::MultiPolygon) = GeoInterface.MultiPolygonTrait()
GeoInterface.geomtrait(::Ngon) = GeoInterface.PolygonTrait()
GeoInterface.geomtrait(::AbstractMesh) = GeoInterface.PolyhedralSurfaceTrait()

GeoInterface.geomtrait(::Simplex{Dim,T,1}) where {Dim,T} = GeoInterface.PointTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,2}) where {Dim,T} = GeoInterface.LineStringTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,3}) where {Dim,T} = GeoInterface.PolygonTrait()

GeoInterface.ncoord(::GeoInterface.PointTrait, g::Point) = length(g)
GeoInterface.getcoord(::GeoInterface.PointTrait, g::Point, i::Int) = g[i]

GeoInterface.ngeom(::GeoInterface.LineTrait, g::Line) = length(g)
GeoInterface.getgeom(::GeoInterface.LineTrait, g::Line, i::Int) = g[i]

GeoInterface.ngeom(::GeoInterface.LineStringTrait, g::LineString) = length(g)
function GeoInterface.getgeom(::GeoInterface.LineStringTrait, g::LineString, i::Int)
    return coordinates(g)[i]
end

GeoInterface.ngeom(::GeoInterface.PolygonTrait, g::Polygon) = length(g.interiors) + 1  # +1 for exterior
function GeoInterface.getgeom(::GeoInterface.PolygonTrait, g::Polygon, i::Int)
    return i > 1 ? g.interiors[i - 1] : g.exterior
end

GeoInterface.ngeom(::GeoInterface.MultiPointTrait, g::MultiPoint) = length(g)
GeoInterface.getgeom(::GeoInterface.MultiPointTrait, g::MultiPoint, i::Int) = g[i]

function GeoInterface.ngeom(::GeoInterface.MultiLineStringTrait, g::MultiLineString)
    return length(g)
end
function GeoInterface.getgeom(::GeoInterface.MultiLineStringTrait, g::MultiLineString,
                              i::Int)
    return g[i]
end

GeoInterface.ngeom(::GeoInterface.MultiPolygonTrait, g::MultiPolygon) = length(g)
GeoInterface.getgeom(::GeoInterface.MultiPolygonTrait, g::MultiPolygon, i::Int) = g[i]

function GeoInterface.ncoord(::GeoInterface.AbstractGeometryTrait,
                             ::Simplex{Dim,T,N,P}) where {Dim,T,N,P}
    return Dim
end
function GeoInterface.ncoord(::GeoInterface.AbstractGeometryTrait,
                             ::AbstractGeometry{Dim,T}) where {Dim,T}
    return Dim
end
function GeoInterface.ngeom(::GeoInterface.AbstractGeometryTrait,
                            ::Simplex{Dim,T,N,P}) where {Dim,T,N,P}
    return N
end
GeoInterface.ngeom(::PolygonTrait, ::Ngon) = 1  # can't have any holes
GeoInterface.getgeom(::PolygonTrait, g::Ngon, _) = LineString(g.points)

GeoInterface.ncoord(::PolyhedralSurfaceTrait, ::Mesh{Dim,T,E,V} where {Dim,T,E,V}) = Dim
GeoInterface.ngeom(::PolyhedralSurfaceTrait, g::AbstractMesh) = length(g)
GeoInterface.getgeom(::PolyhedralSurfaceTrait, g::AbstractMesh, i) = g[i]
