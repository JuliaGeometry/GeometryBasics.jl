"""
    HyperSphere{N, T}

A `HyperSphere` is a generalization of a sphere into N-dimensions.
A `center` and radius, `r`, must be specified.
"""
struct HyperSphere{N,T} <: GeometryPrimitive{N,T}
    center::Point{N,T}
    r::T
end

"""
    Circle{T}

An alias for a HyperSphere of dimension 2. (i.e. `HyperSphere{2, T}`)
"""
const Circle{T} = HyperSphere{2,T}

"""
    Sphere{T}

An alias for a HyperSphere of dimension 3. (i.e. `HyperSphere{3, T}`)
"""
const Sphere{T} = HyperSphere{3,T}

HyperSphere{N}(p::Point{N,T}, number) where {N,T} = HyperSphere{N,T}(p, convert(T, number))

widths(c::HyperSphere{N,T}) where {N,T} = Vec{N,T}(radius(c) * 2)
radius(c::HyperSphere) = c.r
origin(c::HyperSphere) = c.center

Base.minimum(c::HyperSphere{N,T}) where {N,T} = Vec{N,T}(origin(c)) - Vec{N,T}(radius(c))
Base.maximum(c::HyperSphere{N,T}) where {N,T} = Vec{N,T}(origin(c)) + Vec{N,T}(radius(c))

function Base.in(x::Point, c::HyperSphere)
    return norm(origin(c) - x) ≤ radius(c)
end

centered(S::Type{HyperSphere{N,T}}) where {N,T} = S(Vec{N,T}(0), T(0.5))
function centered(::Type{T}) where {T<:HyperSphere}
    return centered(HyperSphere{ndims_or(T, 3),eltype_or(T, Float32)})
end

function coordinates(s::Circle, nvertices=64)
    r = radius(s); o = origin(s)
    ps = [r * Point(cos(phi), sin(phi)) + o for phi in LinRange(0, 2pi, nvertices+1)]
    ps[end] = o
    return ps
end

function texturecoordinates(::Circle, nvertices=64)
    return coordinates(Circle(Point2f(0.5), 0.5f0), nvertices)
end

function faces(::Circle, nvertices=64)
    return [GLTriangleFace(nvertices+1, i, mod1(i+1, nvertices)) for i in 1:nvertices]
end


function coordinates(s::Sphere, nvertices=24)
    θ = LinRange(0, pi, nvertices)
    φ = LinRange(0, 2pi, nvertices)
    inner(θ, φ) = Point(cos(φ) * sin(θ), sin(φ) * sin(θ), cos(θ)) .* s.r .+ s.center
    return [inner(θ, φ) for φ in φ for θ in θ]
end

function texturecoordinates(::Sphere, nvertices=24)
    ux = LinRange(0, 1, nvertices)
    return ivec(((φ, θ) for θ in reverse(ux), φ in ux))
end

function faces(::Sphere, nvertices=24)
    return faces(Rect(0, 0, 1, 1), (nvertices, nvertices))
end

function normals(::Sphere{T}, nvertices=24) where {T}
    return coordinates(Sphere(Point{3,T}(0), 1), nvertices)
end
