# Implementation of trait based interface from https://github.com/JuliaGeo/GeoInterface.jl/

GeoInterface.isgeometry(::Type{<:AbstractGeometry}) = true
GeoInterface.isgeometry(::Type{<:AbstractFace}) = true

GeoInterface.isgeometry(::Type{<:Point}) = true
GeoInterface.isgeometry(::Type{<:AbstractMesh}) = true
GeoInterface.isgeometry(::Type{<:AbstractPolygon}) = true
GeoInterface.isgeometry(::Type{<:LineString}) = true
GeoInterface.isgeometry(::Type{<:MultiPoint}) = true
GeoInterface.isgeometry(::Type{<:MultiLineString}) = true
GeoInterface.isgeometry(::Type{<:MultiPolygon}) = true
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
geointerface_geomtype(::GeoInterface.LineTrait) = Line
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


GeoInterface.ngeom(::LineStringTrait, g::LineString) = length(g)  # n connected points
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
function GeoInterface.getgeom(::MultiLineStringTrait, g::MultiLineString, i::Int)
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
    x, y = GeoInterface.x(geom), GeoInterface.y(geom)
    if GeoInterface.is3d(geom)
        z = GeoInterface.z(geom)
        T = promote_type(typeof(x), typeof(y), typeof(z))
        return Point{3,T}(x, y, z)
    else
        GeoInterface.x(geom), GeoInterface.y(geom)
        T = promote_type(typeof(x), typeof(y))
        return Point{2,T}(x, y)
    end
end

# without a function barrier you get a lot of allocations from runtime types
function _collect_with_type(::Type{PT}, geom) where {PT <: Point{2}}
    return [PT(GeoInterface.x(p), GeoInterface.y(p)) for p in getgeom(geom)]
end

function GeoInterface.convert(::Type{LineString}, type::LineStringTrait, geom)
    g1 = getgeom(geom, 1)
    x, y = GeoInterface.x(g1), GeoInterface.y(g1)
    if GeoInterface.is3d(geom)
        z = GeoInterface.z(g1)
        T = promote_type(typeof(x), typeof(y), typeof(z))
        return LineString([Point{3,T}(GeoInterface.x(p), GeoInterface.y(p), GeoInterface.z(p)) for p in getgeom(geom)])
    else
        T = promote_type(typeof(x), typeof(y))
        return LineString(_collect_with_type(Point{2, T}, geom))
    end
end

function GeoInterface.convert(::Type{Polygon}, type::PolygonTrait, geom)
    t = LineStringTrait()
    exterior = GeoInterface.convert(LineString, t, GeoInterface.getexterior(geom))
    if GeoInterface.nhole(geom) == 0
        return Polygon(exterior)
    else
        interiors = map(h -> GeoInterface.convert(LineString, t, h), GeoInterface.gethole(geom))
        return Polygon(exterior, interiors)
    end
end

function GeoInterface.convert(::Type{MultiPoint}, type::MultiPointTrait, geom)
    g1 = getgeom(geom, 1)
    x, y = GeoInterface.x(g1), GeoInterface.y(g1)
    if GeoInterface.is3d(geom)
        z = GeoInterface.z(g1)
        T = promote_type(typeof(x), typeof(y), typeof(z))
        return MultiPoint([Point{3,T}(GeoInterface.x(p), GeoInterface.y(p), GeoInterface.z(p)) for p in getgeom(geom)])
    else
        T = promote_type(typeof(x), typeof(y))
        return MultiPoint([Point{2,T}(GeoInterface.x(p), GeoInterface.y(p)) for p in getgeom(geom)])
    end
end

function GeoInterface.convert(::Type{MultiLineString}, type::MultiLineStringTrait, geom)
    t = LineStringTrait()
    return MultiLineString(map(l -> GeoInterface.convert(LineString, t, l), getgeom(geom)))
end

function GeoInterface.convert(::Type{MultiPolygon}, type::MultiPolygonTrait, geom)
    t = PolygonTrait()
    return MultiPolygon(map(poly -> GeoInterface.convert(Polygon, t, poly), getgeom(geom)))
end

function Extents.extent(rect::Rect2)
    (xmin, ymin), (xmax, ymax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end

function Extents.extent(rect::Rect3)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax), Z = (zmin, zmax))
end

function Extents.extent(rect::Rect4)
    (xmin, ymin, zmin, mmin), (xmax, ymax, zmax, mmax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax), Z = (zmin, zmax), M = (mmin, mmax))
end
