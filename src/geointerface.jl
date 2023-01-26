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

GeoInterface.geomtrait(::Point) = PointTrait()
GeoInterface.geomtrait(::Line) = LineTrait()
GeoInterface.geomtrait(::LineString) = LineStringTrait()
GeoInterface.geomtrait(::Polygon) = PolygonTrait()
GeoInterface.geomtrait(::MultiPoint) = MultiPointTrait()
GeoInterface.geomtrait(::MultiLineString) = MultiLineStringTrait()
GeoInterface.geomtrait(::MultiPolygon) = MultiPolygonTrait()
GeoInterface.geomtrait(::Ngon) = PolygonTrait()
GeoInterface.geomtrait(::AbstractMesh) = PolyhedralSurfaceTrait()

GeoInterface.geomtrait(::Simplex{Dim,T,1}) where {Dim,T} = PointTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,2}) where {Dim,T} = LineStringTrait()
GeoInterface.geomtrait(::Simplex{Dim,T,3}) where {Dim,T} = PolygonTrait()

# GeoInterface calls this method in `GeoInterface.convert(GeometryBasics, ...)`
geointerface_geomtype(::GeoInterface.PointTrait) = Point
geointerface_geomtype(::GeoInterface.MultiPointTrait) = MultiPoint
geointerface_geomtype(::GeoInterface.LineStringTrait) = LineString
geointerface_geomtype(::GeoInterface.MultiLineStringTrait) = MultiLineString
geointerface_geomtype(::GeoInterface.PolygonTrait) = Polygon
geointerface_geomtype(::GeoInterface.MultiPolygonTrait) = MultiPolygon
geointerface_geomtype(::GeoInterface.PolyhedralSurfaceTrait) = AbstractMesh

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
                              i::Int)::typeof(g.exterior)
    return i > 1 ? g.interiors[i - 1] : g.exterior
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
                             ::Simplex{Dim,T,N,P}) where {Dim,T,N,P}
    return Dim
end
function GeoInterface.ncoord(::AbstractGeometryTrait,
                             ::AbstractGeometry{Dim,T}) where {Dim,T}
    return Dim
end
function GeoInterface.ngeom(::AbstractGeometryTrait,
                            ::Simplex{Dim,T,N,P}) where {Dim,T,N,P}
    return N
end
GeoInterface.ngeom(::PolygonTrait, ::Ngon) = 1  # can't have any holes
GeoInterface.getgeom(::PolygonTrait, g::Ngon, _) = LineString(g.points)

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

GeoInterface.centroid(::Union{GeoInterface.MultiPolygonTrait, GeoInterface.PolygonTrait}, geom::Union{AbstractGeometry, MultiPolygon, MultiLineString}) = centroid(geom)


# Implementations for and overloads of various GeoInterface optional functions

function GeoInterface.extent(::GeoInterface.AbstractGeometryTrait, a::Union{GeometryBasics.AbstractGeometry, GeometryBasics.GeometryPrimitive, AbstractVector{<: GeometryBasics.AbstractGeometry}})
    bbox = Rect(a)
    return if ndims(bbox) == 1
        GeoInterface.Extents.Extent(X = (minimum(bbox), maximum(bbox)))
    elseif ndims(bbox) == 2
        xmin, ymin = minimum(bbox)
        xmax, ymax = maximum(bbox)
        GeoInterface.Extents.Extent(X = (xmin, xmax), Y = (ymin, ymax))
    else # ndims ≥ 3
        xmin, ymin, zmin = minimum(bbox)
        xmax, ymax, zmax = maximum(bbox)
        GeoInterface.Extents.Extent(X = (xmin, xmax), Y = (ymin, ymax), Z = (zmin, zmax))
    end
end

# quick and dirty conversions from GeoInterface's returned Extent back to Rect
function Rect{2, T}(ext::GeoInterface.Extents.Extent{(:X, :Y)}) where {T}
    xmin, xmax = ext.X
    ymin, ymax = ext.Y
    return Rect{2, T}(xmin, ymin, xmax - xmin, ymax - ymin)
end

function Rect{3, T}(ext::GeoInterface.Extents.Extent{(:X, :Y, :Z)}) where {T}
    xmin, xmax = ext.X
    ymin, ymax = ext.Y
    zmin, zmax = ext.Z
    return Rect{3, T}(xmin, ymin, zmin, xmax - xmin, ymax - ymin, zmax - zmin)
end

function signed_area(a::Point{2, T}, b::Point{2, T}, c::Point{2, T}) where T
    return ((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
end

function signed_area(points::AbstractVector{<: Point{2, T}}) where {T}
    area = sum((points[i][1] * points[i+1][2] - points[i][2] * points[i+1][1] for i in 1:(length(points)-1))) / 2.0
end

function signed_area(ls::GeometryBasics.LineString)
    # coords = GeometryBasics.decompose(Point2f, ls)
    return sum((p1[1] * p2[2] - p1[2] * p2[1] for (p1, p2) in ls)) / 2.0#signed_area(coords)
end

function signed_area(poly::GeometryBasics.Polygon{2})
    area = abs(signed_area(poly.exterior))
    for hole in poly.interiors
        area -= abs(signed_area(hole))
    end
    return area
end

signed_area(mp::MultiPolygon) = sum(signed_area.(mp.polygons))

function centroid(ls::LineString{2, T}) where T
    centroid = Point{2, T}(0)
    total_area = T(0)
    if length(ls) == 1
        return sum(ls[1])/2
    end

    p0 = ls[1][1]

    for i in 1:(length(ls)-1)
        p1 = ls[i][2]
        p2 = ls[i+1][2]
        area = signed_area(p0, p1, p2)
        centroid = centroid .+ Point{2, T}((p0[1] + p1[1] + p2[1])/3, (p0[2] + p1[2] + p2[2])/3) * area
        total_area += area
    end
    return centroid ./ total_area
end

function centroid(poly::GeometryBasics.Polygon{2, T}) where T
    exterior_points = decompose(Point2f, poly.exterior)
    exterior_centroid = centroid(poly.exterior)
    exterior_area = signed_area(poly.exterior)

    total_area = exterior_area
    interior_numerator = Point{2, T}(0)
    for interior in poly.interiors
        interior_points = decompose(Point2f, interior)
        interior_centroid = centroid(interior)
        interior_area = signed_area(interior_points)
        total_area += interior_area
        interior_numerator += interior_centroid * interior_area
    end

    return (exterior_centroid * exterior_area - interior_numerator) / total_area
        
end

function centroid(multipoly::MultiPolygon)

    centroids = centroid.(multipoly.polygons)

    areas = signed_area.(multipoly.polygons)
    areas ./= sum(areas)

    return sum(centroids .* areas) / sum(areas)

end


function centroid(rect::Rect{N, T}) where {N, T}
    return Point{N, T}(rect.origin .- rect.widths ./ 2)
end

# function distance(poly::Polygon{N, T1}, point::Point{N, T2}) where {N, T1, T2}
#     FinalType = promote_type(T1, T2)
# end

function contains(ls::GeometryBasics.LineString{2, T1}, point::Point{2, T2}) where {T1, T2}

    # the original C code from https://wrfranklin.org/Research/Short_Notes/pnpoly.html
    # int pnpoly(int npol, float *xp, float *yp, float x, float y)
    # {
    #   i = 1
    #   j = 1
    #   c = false
    #   for (i = 0, j = npol-1; i < npol; j = i++) {
    #     if ((((yp[i]<=y) && (y<yp[j])) ||
    #          ((yp[j]<=y) && (y<yp[i]))) &&
    #         (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))

    #       c = !c;
    #   }
    #   return c;

    # }

    x, y = point

    c = false
    for (p1, p2) in ls
        if ((p1[2] ≤ y && y < p2[2]) ||
            (p2[2] ≤ y && y < p1[2])) &&
            (x < (p2[1] - p1[1]) * (y - p1[2]) / (p2[2] - p1[2]) + p1[1])
            c = !c
        end
    end

    return c

end

function contains(poly::Polygon{2, T1}, point::Point{2, T2}) where {T1, T2}
    c = contains(poly.exterior, point)
    for interior in poly.interiors
        if contains(interior, point)
            c = false # if hole contains point, then point is not contained in poly.
            break
        end
    end
    return c
end

contains(mp::MultiPolygon{2, T1}, point::Point{2, T2}) where {T1, T2} = any((contains(poly, point) for poly in polygons))

GeoInterface.contains(::GeoInterface.LineStringTrait, ::GeoInterface.PointTrait, ls::LineString{2, T1}, point::Point{2, T2}) where {T1, T2} = contains(ls, point)
GeoInterface.contains(::GeoInterface.PolygonTrait, ::GeoInterface.PointTrait, poly::Polygon{2, T1}, point::Point{2, T2}) where {T1, T2} = contains(poly, point)
GeoInterface.contains(::GeoInterface.MultiPolygonTrait, ::GeoInterface.PointTrait, multipoly::MultiPolygon{2, T1}, point::Point{2, T2}) where {T1, T2} = contains(multipoly, point)