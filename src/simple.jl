
"""
A `HyperRectangle` is a generalization of a rectangle into N-dimensions.
Formally it is the cartesian product of intervals, which is represented by the
`origin` and `width` fields, whose indices correspond to each of the `N` axes.
"""
struct HyperRectangle{N, T} <: GeometryPrimitive{N, T}
    origin::Vec{N, T}
    widths::Vec{N, T}
end


"""
A `HyperSphere` is a generalization of a sphere into N-dimensions.
A `center` and radius, `r`, must be specified.
"""
struct HyperSphere{N, T} <: GeometryPrimitive{N, T}
    center::Point{N, T}
    r::T
end
const Circle{T} = HyperSphere{2, T}
const Sphere{T} = HyperSphere{3, T}

widths(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(radius(c) * T(2))
radius(c::HyperSphere) = c.r
origin(c::HyperSphere) = c.center
minimum(c::HyperSphere) = origin(c) .- radius(c)
maximum(c::HyperSphere) = origin(c) .+ radius(c)

function Base.in(point::Point{N, T}, c::HyperSphere{N}) where {N, T}
    xΔ = norm(origin(c) .- point)
    return xΔ <= radius(c)
end

centered(S::Type{HyperSphere{N, T}}) where {N, T} = S(Vec{N, T}(0), T(0.5))
centered(::Type{T}) where {T <: HyperSphere} = centered(HyperSphere{ndims_or(T, 3), eltype_or(T, Float32)})


struct Pyramid{T} <: GeometryPrimitive{3, T}
    middle::Point{3, T}
    length::T
    width ::T
end

"""
A `Cylinder` is a 2D rectangle or a 3D cylinder defined by its origin point,
its extremity and a radius. `origin`, `extremity` and `r`, must be specified.
"""
struct Cylinder{N,T<: AbstractFloat} <: GeometryPrimitive{N,T}
    origin::Point{N,T}
    extremity::Point{N,T}
    r::T
end

origin(c::Cylinder) = c.origin
extremity(c::Cylinder) = c.extremity
radius(c::Cylinder) = c.r
height(c::Cylinder) = norm(c.extremity - c.origin)
direction(c::Cylinder) = (c.extremity .- c.origin) ./ height(c)

function rotation(c::Cylinder{2, T}) where T
    d2 = direction(c); u = Vec(d2[1], d2[2], T(0))
    v = normalize(Vec(u[2], -u[1], T(0)))
    return hcat(v, u, Vec{3, T}(0, 0, 1))
end

function rotation(c::Cylinder{3, T}) where T
    d3 = direction(c); u = Vec(d3[1], d3[2], d3[3])
    if abs(u[1]) > 0 || abs(u[2]) > 0
        v = Vec(u[2], -u[1], T(0))
    else
        v = Vec(T(0), -u[3], u[2])
    end
    v = normalize(v)
    w = Vec(u[2] * v[3] - u[3] * v[2], -u[1] * v[3] + u[3] * v[1], u[1] * v[2] - u[2] * v[1])
    return hcat(v, w, u)
end
