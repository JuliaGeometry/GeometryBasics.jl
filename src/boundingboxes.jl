"""
    boundbox(geometry)

Axis-aligned bounding box of the `geometry`.
"""
boundbox(geom) = boundbox(coordinates(geom))

# fallback implementation treats geometry as
# a set of points (i.e. coordinates)
function boundbox(geometry::AbstractArray{<:AbstractPoint{N,T}}) where {N,T}
    vmin = Point{N,T}(typemax(T))
    vmax = Point{N,T}(typemin(T))
    for p in geometry
        vmin, vmax = minmax(p, vmin, vmax)
    end
    Rect{N,T}(vmin, vmax - vmin)
end

# --------------
# SPECIAL CASES
# --------------

boundbox(a::Rect) = a

function boundbox(a::Pyramid{T}) where {T}
    w = a.width / 2
    h = a.length
    m = a.middle
    Rect{3,T}(m - Point{3,T}(w, w, 0), m + Point{3,T}(w, w, h))
end

function boundbox(a::Sphere{T}) where {T}
    mini, maxi = extrema(a)
    Rect{3,T}(mini, maxi .- mini)
end
