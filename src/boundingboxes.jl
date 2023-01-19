# the basic method
"""
Construct a HyperRectangle enclosing all points.
"""
function Rect{N1,T1}(geometry::AbstractArray{PT}) where {N1,T1,PT<:AbstractPoint}
    N2, T2 = length(PT), eltype(PT)
    @assert N1 >= N2
    vmin = Point{N2,T2}(typemax(T2))
    vmax = Point{N2,T2}(typemin(T2))
    for p in geometry
        vmin, vmax = minmax(p, vmin, vmax)
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

function Rect(geometry::AbstractArray{<:Point{N,T}}) where {N, T}
    return Rect{N,T}(geometry)
end

function Rect{T1}(geometry::AbstractArray{<:Point{N,T2}}) where {N, T1, T2}
    return Rect{N,T1}(geometry)
end

Rect(a::AbstractGeometry{N, T}) where {N, T} = Rect{N, T}(a)
Rect{T1}(a::AbstractGeometry{N, T2}) where {N, T1, T2} = Rect{N, T1}(a)
function Rect{N1, T1}(a::AbstractGeometry{N2, T2}) where {N1, T1, N2, T2}
    return Rect{N1, T1}(coordinates(a))
end

Rect(a::AbstractArray{<: AbstractGeometry{N, T}}) where {N, T} = Rect{N, T}(a)
Rect{T1}(a::AbstractArray{<: AbstractGeometry{N, T2}}) where {N, T1, T2} = Rect{N, T1}(a)
function Rect{N1, T1}(a::AbstractArray{<: AbstractGeometry{N2, T2}}) where {N1, T1, N2, T2}
    return reduce(union, Rect{N1, T1}.(a))
end

# generic fallbacks and entry point
Rect{T}(a) where {T} = Rect{T}(coordinates(a))
Rect{N,T}(a) where {N,T} = Rect{N,T}(coordinates(a))

function Rect(primitive::GeometryPrimitive{N,T}) where {N,T}
    return Rect{N,T}(primitive)
end

function Rect{T}(primitive::GeometryPrimitive{N,T}) where {N,T}
    return Rect{N,T}(primitive)
end

# function Rect{N1, T1}(primitives::AbstractArray{<: GeometryPrimitive{N2, T2}})
# end

# specialized overloads for geometries with a known bbox

function Rect{3, T}(a::Pyramid) where {T}
    w, h = a.width / T(2), a.length
    m = Vec{3,T}(a.middle)
    return Rect{3, T}(m .- Vec{3,T}(w, w, 0), m .+ Vec{3,T}(w, w, h))
end

function Rect{2, T}(a::Sphere) where {T}
    mini, maxi = extrema(a)
    return Rect{T}(mini, maxi .- mini)
end

function Rect{N1, T}(hypersph::GeometryBasics.HyperSphere{N2}) where {T, N1, N2}
    mini, maxi = extrema(hypersph)
    if length(mini) < N1
        mini = vcat(mini, zeros(N1 - length(mini)))
        maxi = vcat(maxi, zeros(N1 - length(mini)))
    end
    return Rect{N1, T}(mini, maxi .- mini)
end
