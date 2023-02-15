# Implementation of trait based interface from https://github.com/JuliaGeo/GeoInterface.jl/

GeoInterface.isgeometry(::Type{<:AbstractGeometry}) = true
GeoInterface.isgeometry(::Type{<:AbstractFace}) = true
GeoInterface.isgeometry(::Type{<:Point}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractGeometry}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:Point}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:LineString}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractPolygon}}) = true
GeoInterface.isgeometry(::Type{<:AbstractVector{<:AbstractFace}}) = true
GeoInterface.isgeometry(::Type{<:Mesh}) = true

GeoInterface.geomtrait(::Point) = PointTrait()
GeoInterface.geomtrait(::Line) = LineTrait()
GeoInterface.geomtrait(::LineString) = LineStringTrait()
GeoInterface.geomtrait(::Polygon) = PolygonTrait()
GeoInterface.geomtrait(::MultiPoint) = MultiPointTrait()
GeoInterface.geomtrait(::MultiLineString) = MultiLineStringTrait()
GeoInterface.geomtrait(::MultiPolygon) = MultiPolygonTrait()
GeoInterface.geomtrait(::Ngon) = PolygonTrait()
GeoInterface.geomtrait(::AbstractMesh) = PolyhedralSurfaceTrait()

# GeoInterface calls this method in `GeoInterface.convert(GeometryBasics, ...)`
geointerface_geomtype(::GeoInterface.PointTrait) = Point
geointerface_geomtype(::GeoInterface.MultiPointTrait) = MultiPoint
geointerface_geomtype(::GeoInterface.LineStringTrait) = LineString
geointerface_geomtype(::GeoInterface.MultiLineStringTrait) = MultiLineString
geointerface_geomtype(::GeoInterface.PolygonTrait) = Polygon
geointerface_geomtype(::GeoInterface.MultiPolygonTrait) = MultiPolygon
geointerface_geomtype(::GeoInterface.PolyhedralSurfaceTrait) = Mesh

GeoInterface.geomtrait(::Simplex{Dim,T,1}) where {Dim,T} = PointTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,2}) where {Dim,T} = LineStringTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,3}) where {Dim,T} = PolygonTrait()

GeoInterface.ncoord(::PointTrait, g::Point) = length(g)
GeoInterface.getcoord(::PointTrait, g::Point, i::Int) = g[i]

GeoInterface.ngeom(::LineTrait, g::Line) = length(g)
GeoInterface.getgeom(::LineTrait, g::Line, i::Int) = g[i]


GeoInterface.ngeom(::LineStringTrait, g::LineString) = length(g) + 1  # n line segments + 1
GeoInterface.ncoord(::LineStringTrait, g::LineString{Dim}) where {Dim} = Dim
function GeoInterface.getgeom(::LineStringTrait, g::LineString, i::Int)
    return GeometryBasics.coordinates(g)[i]
end

GeoInterface.ngeom(::PolygonTrait, g::Polygon) = length(g.interiors) + 1  # +1 for exterior
function GeoInterface.getgeom(::PolygonTrait,
                              g::Polygon,
                              i::Int)
    return i > 1 ? LineString(g.interiors[i - 1]) : LineString(g.exterior)
end

GeoInterface.ngeom(::MultiPointTrait, g::MultiPoint) = length(g)
GeoInterface.getgeom(::MultiPointTrait, g::MultiPoint, i::Int) = g[i]

function GeoInterface.ngeom(::MultiLineStringTrait, g::MultiLineString)
    return length(g)
end
function GeoInterface.getgeom(::MultiLineStringTrait, g::MultiLineString,
                              i::Int)
    return g[i]
end
GeoInterface.ncoord(::MultiLineStringTrait, g::MultiLineString{Dim}) where {Dim} = Dim

GeoInterface.ngeom(::MultiPolygonTrait, g::MultiPolygon) = length(g)
GeoInterface.getgeom(::MultiPolygonTrait, g::MultiPolygon, i::Int) = g[i]

function GeoInterface.ncoord(::AbstractGeometryTrait,
                             ::Simplex{Dim,T,N}) where {Dim,T,N}
    return Dim
end
function GeoInterface.ncoord(::AbstractGeometryTrait,
                             ::AbstractGeometry{Dim,T}) where {Dim,T}
    return Dim
end
function GeoInterface.ngeom(::AbstractGeometryTrait,
                            ::Simplex{Dim,T,N}) where {Dim,T,N}
    return N
end
GeoInterface.ngeom(::PolygonTrait, ::Ngon) = 1  # can't have any holes
GeoInterface.getgeom(::PolygonTrait, g::Ngon, _) = LineString([g.points...])

function GeoInterface.ncoord(::PolyhedralSurfaceTrait,
                             ::Mesh{Dim,T,E,V} where {Dim,T,E,V})
    return Dim
end
GeoInterface.ngeom(::PolyhedralSurfaceTrait, g::AbstractMesh) = length(g)
GeoInterface.getgeom(::PolyhedralSurfaceTrait, g::AbstractMesh, i) = g[i]

function GeoInterface.convert(::Type{Point}, type::PointTrait, geom)
    dim = Int(ncoord(geom))
    return Point{dim, Float64}(GeoInterface.coordinates(geom))
end

function GeoInterface.convert(::Type{LineString}, type::LineStringTrait, geom)
    dim = Int(ncoord(geom))
    return LineString([Point{dim, Float64}(GeoInterface.coordinates(p)) for p in getgeom(geom)])
end

function GeoInterface.convert(::Type{Polygon}, type::PolygonTrait, geom)
    t = LineStringTrait()
    exterior = GeoInterface.convert(LineString, t, GeoInterface.getexterior(geom))
    if GeoInterface.nhole(geom) == 0
        return Polygon(exterior)
    else
        interiors = GeoInterface.convert.(LineString, Ref(t), GeoInterface.gethole(geom))
        return Polygon(exterior, interiors)
    end
end

function GeoInterface.convert(::Type{MultiPoint}, type::MultiPointTrait, geom)
    dim = Int(ncoord(geom))
    return MultiPoint([Point{dim, Float64}(GeoInterface.coordinates(p)) for p in getgeom(geom)])
end

function GeoInterface.convert(::Type{MultiLineString}, type::MultiLineStringTrait, geom)
    t = LineStringTrait()
    return MultiLineString([GeoInterface.convert(LineString, t, l) for l in getgeom(geom)])
end

function GeoInterface.convert(::Type{MultiPolygon}, type::MultiPolygonTrait, geom)
    t = PolygonTrait()
    return MultiPolygon([GeoInterface.convert(Polygon, t, poly) for poly in getgeom(geom)])
end
