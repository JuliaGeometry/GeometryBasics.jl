function Rect(geometry::AbstractArray{<:Point{N,T}}) where {N,T}
    return Rect{N,T}(geometry)
end

"""
    Rect(points::AbstractArray{<: Point})

Construct a bounding box countaining all the given points.
"""
function Rect{N1,T1}(geometry::AbstractArray{PT}) where {N1,T1,PT<:Point}
    N2, T2 = length(PT), eltype(PT)
    @assert N1 >= N2
    vmin = Point{N2,T2}(typemax(T2))
    vmax = Point{N2,T2}(typemin(T2))
    for p in geometry
        vmin, vmax = _minmax(p, vmin, vmax)
    end
    o = vmin
    w = vmax - vmin
    return if N1 > N2
        z = zero(Vec{N1 - N2,T1})
        Rect{N1,T1}(vcat(o, z), vcat(w, z))
    else
        Rect{N1,T1}(o, w)
    end
end

"""
    Rect(primitive::GeometryPrimitive)

Construct a bounding box for the given primitive.
"""
function Rect(primitive::GeometryPrimitive{N,T}) where {N,T}
    return Rect{N,T}(primitive)
end

function Rect{T}(primitive::GeometryPrimitive{N,T}) where {N,T}
    return Rect{N,T}(primitive)
end

function Rect{T}(a::Pyramid) where {T}
    w, h = a.width / T(2), a.length
    m = Vec{3,T}(a.middle)
    return Rect{T}(m .- Vec{3,T}(w, w, 0), m .+ Vec{3,T}(w, w, h))
end

function Rect{T}(a::Sphere) where {T}
    mini, maxi = extrema(a)
    return Rect{T}(mini, maxi .- mini)
end

Rect{T}(a) where {T} = Rect{T}(coordinates(a))
Rect{N,T}(a) where {N,T} = Rect{N,T}(coordinates(a))
