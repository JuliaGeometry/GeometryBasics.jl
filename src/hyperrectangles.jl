origin(prim::HyperRectangle) = prim.origin
maximum(prim::HyperRectangle) = origin(prim) + widths(prim)
minimum(prim::HyperRectangle) = origin(prim)
length(prim::HyperRectangle{N}) where N = N
widths(prim::HyperRectangle) = prim.widths

"""
Splits an HyperRectangle into two along an axis at a given location.
"""
split(b::HyperRectangle, axis, value::Integer) = _split(b, axis, value)
split(b::HyperRectangle, axis, value::Number) = _split(b, axis, value)
function _split(b::H, axis, value) where H<:HyperRectangle
    bmin = minimum(b)
    bmax = maximum(b)
    b1max = setindex(bmax, value, axis)
    b2min = setindex(bmin, value, axis)
    return (H(bmin, b1max-bmin), H(b2min, bmax-b2min))
end

# empty constructor such that update will always include the first point
function Rect{N, T}() where {T, N}
    Rect(Vec{N,T}(typemax(T)), Vec{N,T}(typemin(T)))
end

# conversion from other HyperRectangles
function Rect{N,T1}(a::Rect{N, T2}) where {N,T1,T2}
    Rect(Vec{N, T1}(minimum(a)), Vec{N, T1}(widths(a)))
end

function Rect(v1::Vec{N, T1}, v2::Vec{N, T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    Rect{N,T}(Vec{N, T}(v1), Vec{N, T}(v2))
end
function Rect{N}(v1::Vec{N, T1}, v2::Vec{N, T2}) where {N,T1,T2}
    T = promote_type(T1, T2)
    Rect{N,T}(Vec{N, T}(v1), Vec{N, T}(v2))
end

function Rect{N, T}(a::GeometryPrimitive) where {N, T}
    Rect{N, T}(Vec{N, T}(minimum(a)), Vec{N, T}(widths(a)))
end

"""
```
HyperRectangle(vals::Number...)
```
HyperRectangle constructor for indidually specified intervals.
e.g. HyperRectangle(0,0,1,2) has origin == Vec(0,0) and
width == Vec(1,2)
"""
function Rect(vals::Vararg{Number, N}) where N
    iseven(N) || error("Construct needs even number of arguments! Found: $vals")
    mid = div(N, 2)
    origin = ntuple(i-> vals[i], mid)
    widths = ntuple(i-> vals[mid + i], mid)
    Rect{mid}(origin, widths)
end

#=
From other types
=#
function Rect{2}(xy::NamedTuple{(:x, :y)}, wh::NamedTuple{(:width, :height)})
    Rect{2}(values(xy), values(wh))
end

function Rect{3}(x::Rect{2, T}) where T
    Rect{3, T}(Vec{3, T}(minimum(x)..., 0), Vec{3, T}(widths(x)..., 0.0))
end

function Rect{2}(xy::VecTypes{2}, w::Number, h::Number)
    Rect{2}(xy..., w, h)
end

function Rect{2}(x::Number, y::Number, wh::VecTypes{2})
    Rect{2}(x, y, wh...)
end

#=
From limits
=#
function Rect{3}(x::Tuple{Tuple{<: Number, <: Number}, Tuple{<: Number, <: Number}})
    Rect{3}(Vec3f0(x[1]..., 0), Vec3f0(x[2]..., 0))
end
function Rect{3}(x::Tuple{Tuple{<: Number, <: Number, <: Number}, Tuple{<: Number, <: Number, <: Number}})
    Rect{3}(Vec3f0(x[1]...), Vec3f0(x[2]...))
end


"""
Transform a `HyperRectangle` using a matrix. Maintains axis-align properties
so a significantly larger HyperRectangle may be generated.
"""
function *(m::Mat{N1, N1, T1}, h::Rect{N2, T2}) where {N1, N2, T1, T2}
    # TypeVar constants
    T = promote_type(T1, T2)
    D = N1 - N2

    # get all points on the HyperRectangle
    d = decompose(Point, h)
    # make sure our points are sized for the tranform
    pts = ntuple(i-> Vec{N1, T}(vcat(d[i], ones(Vec{D, T}))), 2^N2)::NTuple{2^N2,Vec{N1,T}}

    vmin = Vec{N1, T}(typemax(T))
    vmax = Vec{N1, T}(typemin(T))
    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * pt
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    Rect(vmin, vmax .- vmin)
end


function -(h::Rect{N}, move::Vec{N}) where N
    Rect(minimum(h) .- move, widths(h))
end


function *(m::Mat{N, N, T1}, h::Rect{N, T2}) where {N,T1,T2}
    # equal dimension case

    # TypeVar constants
    T = promote_type(T1, T2)

    # get all points on the HyperRectangle
    pts = decompose(Point, h)

    # make sure our points are sized for the tranform
    vmin = Vec{N, T}(typemax(T))
    vmax = Vec{N, T}(typemin(T))

    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * Vec(pt)
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    HyperRectangle{N,T}(vmin, vmax-vmin)
end

# fast path. TODO make other versions fast without code duplications like now
function *(m::Mat{4,4,T}, h::Rect{3,T}) where T
    # equal dimension case

    # get all points on the HyperRectangle
    pts = (
        Vec{4,T}(0.0,0.0,0.0, 1.0),
        Vec{4,T}(1.0,0.0,0.0, 1.0),
        Vec{4,T}(0.0,1.0,0.0, 1.0),
        Vec{4,T}(1.0,1.0,0.0, 1.0),
        Vec{4,T}(0.0,0.0,1.0, 1.0),
        Vec{4,T}(1.0,0.0,1.0, 1.0),
        Vec{4,T}(0.0,1.0,1.0, 1.0),
        Vec{4,T}(1.0,1.0,1.0, 1.0)
    )

    # make sure our points are sized for the tranform
    vmin = Vec{4, T}(typemax(T))
    vmax = Vec{4, T}(typemin(T))
    o, w = origin(h), widths(h)
    _o = Vec{4, T}(o[1], o[2], o[3], T(0))
    _w = Vec{4, T}(w[1], w[2], w[3], T(1))
    # tranform all points, tracking min and max points
    for pt in pts
        pn = m * (_o + (pt .* _w))
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    _vmin = Vec{3, T}(vmin[1], vmin[2], vmin[3])
    _vmax = Vec{3, T}(vmax[1], vmax[2], vmax[3])
    Rect(_vmin, _vmax - _vmin)
end

function Rect(geometry::AbstractArray{<: Point{N, T}}) where {N,T}
    Rect{N,T}(geometry)
end

@inline function minmax(p::StaticVector, vmin, vmax)
    isnan(p) && return (vmin, vmax)
    min.(p, vmin), max.(p, vmax)
end

# Annoying special case for view(Vector{Point}, Vector{Face})
@inline function minmax(tup::Tuple, vmin, vmax)
    for p in tup
        isnan(p) && continue
        vmin = min.(p, vmin)
        vmax = max.(p, vmax)
    end
    vmin, vmax
end

"""
Construct a HyperRectangle enclosing all points.
"""
function Rect{N1, T1}(
        geometry::AbstractArray{PT}
    ) where {N1, T1, PT <: Point}
    N2, T2 = length(PT), eltype(PT)
    @assert N1 >= N2
    vmin = Point{N2, T2}(typemax(T2))
    vmax = Point{N2, T2}(typemin(T2))
    for p in geometry
        vmin, vmax = minmax(p, vmin, vmax)
    end
    o = vmin; w = vmax - vmin
    if N1 > N2
        z = zero(Vec{N1-N2, T1})
        return Rect(vcat(o, z), vcat(w, z))
    else
        return Rect(o, w)
   end
end

function Rect{T}(a::Pyramid) where T
    w,h = a.width/T(2), a.length
    m = Vec{3,T}(a.middle)
    Rect{T}(m-Vec{3,T}(w,w,0), m+Vec{3,T}(w, w, h))
end

Rect{T}(a::Cube) where T = Rect{T}(origin(a), widths(a))

Rect{T}(a::AbstractMesh) where T = Rect{T}(vertices(a))

function positive_widths(rect::Rect{N, T}) where {N, T}
    mini, maxi = minimum(rect), maximum(rect)
    realmin = min.(mini, maxi)
    realmax = max.(mini, maxi)
    Rect{N, T}(realmin, realmax .- realmin)
end


# set operations

"""
Perform a union between two HyperRectangles.
"""
function union(h1::Rect{N}, h2::Rect{N}) where N
    m = min.(minimum(h1), minimum(h2))
    mm = max.(maximum(h1), maximum(h2))
    Rect{N}(m, mm - m)
end


"""
Perform a difference between two HyperRectangles.
"""
diff(h1::Rect, h2::Rect) = h1


"""
Perform a intersection between two HyperRectangles.
"""
function intersect(h1::Rect{N}, h2::Rect{N}) where N
    m = max.(minimum(h1), minimum(h2))
    mm = min.(maximum(h1), maximum(h2))
    Rect{N}(m, mm - m)
end

function intersect(a::SimpleRectangle, b::SimpleRectangle)
    min_n = max.(minimum(a), minimum(b))
    max_n = min.(maximum(a), maximum(b))
    w = max_n - min_n
    SimpleRectangle(min_n[1], min_n[2], w[1], w[2])
end



function update(b::HyperRectangle{N, T}, v::Vec{N, T2}) where {N, T, T2}
    update(b, Vec{N, T}(v))
end
function update(b::HyperRectangle{N, T}, v::Vec{N, T}) where {N, T}
    m = min.(minimum(b), v)
    maxi = maximum(b)
    mm = if isnan(maxi)
        v-m
    else
        max.(v, maxi) - m
    end
    HyperRectangle{N, T}(m, mm)
end

# Min maximum distance functions between hrectangle and point for a given dimension
@inline function min_dist_dim(rect::HyperRectangle{N, T}, p::Vec{N, T}, dim::Int) where {N, T}
    max(zero(T), max(minimum(rect)[dim] - p[dim], p[dim] - maximum(rect)[dim]))
end

@inline function max_dist_dim(rect::HyperRectangle{N, T}, p::Vec{N, T}, dim::Int) where {N, T}
    max(maximum(rect)[dim] - p[dim], p[dim] - minimum(rect)[dim])
end

@inline function min_dist_dim(rect1::HyperRectangle{N, T},
                              rect2::HyperRectangle{N, T},
                              dim::Int) where {N, T}
    max(zero(T), max(
        minimum(rect1)[dim] - maximum(rect2)[dim],
        minimum(rect2)[dim] - maximum(rect1)[dim]
    ))
end

@inline function max_dist_dim(rect1::HyperRectangle{N, T},
                              rect2::HyperRectangle{N, T},
                              dim::Int) where {N, T}
    max(
        maximum(rect1)[dim] - minimum(rect2)[dim],
        maximum(rect2)[dim] - minimum(rect1)[dim]
    )
end

# Total minimum maximum distance functions
@inline function min_euclideansq(rect::HyperRectangle{N, T},
                                 p::Union{Vec{N, T},
                                 HyperRectangle{N, T}}) where {N, T}
    minimum_dist = T(0.0)
    for dim in 1:length(p)
        d = min_dist_dim(rect, p, dim)
        minimum_dist += d*d
    end
    return minimum_dist
end


@inline function max_euclideansq(rect::HyperRectangle{N, T}, p::Union{Vec{N, T}, HyperRectangle{N, T}}) where {N, T}
    maximum_dist = T(0.0)
    for dim in 1:length(p)
        d = max_dist_dim(rect, p, dim)
        maximum_dist += d*d
    end
    return maximum_dist
end
function max_euclidean(rect::HyperRectangle{N, T}, p::Union{Vec{N, T}, HyperRectangle{N, T}}) where {N, T}
    sqrt(max_euclideansq(rect, p))
end


# Functions that return both minimum and maximum for convenience
@inline function minmax_dist_dim(rect::HyperRectangle{N, T}, p::Union{Vec{N, T}, HyperRectangle{N, T}}, dim::Int) where {N, T}
    minimum_d = min_dist_dim(rect, p, dim)
    maximum_d = max_dist_dim(rect, p, dim)
    return minimum_d, maximum_d
end


@inline function minmax_euclideansq(rect::HyperRectangle{N, T}, p::Union{Vec{N, T}, HyperRectangle{N, T}}) where {N, T}
    minimum_dist = min_euclideansq(rect, p)
    maximum_dist = max_euclideansq(rect, p)
    return minimum_dist, maximum_dist
end

function minmax_euclidean(rect::HyperRectangle{N, T}, p::Union{Vec{N, T}, HyperRectangle{N, T}}) where {N, T}
    minimumsq, maximumsq = minmax_euclideansq(rect, p)
    return sqrt(minimumsq), sqrt(maximumsq)
end

#
# http://en.wikipedia.org/wiki/Allen%27s_interval_algebra
#

function before(b1::Rect{N}, b2::Rect{N}) where N
    for i = 1:N
        maximum(b1)[i] < minimum(b2)[i] || return false
    end
    true
end

meets(b1::Rect{N}, b2::Rect{N}) where N = maximum(b1) == minimum(b2)

function overlaps(b1::Rect{N}, b2::Rect{N}) where N
    for i = 1:N
        maximum(b2)[i] > maximum(b1)[i] > minimum(b2)[i] &&
        minimum(b1)[i] < minimum(b2)[i] || return false
    end
    true
end

function starts(b1::Rect{N}, b2::Rect{N}) where N
    if minimum(b1) == minimum(b2)
        for i = 1:N
            maximum(b1)[i] < maximum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

function during(b1::Rect{N}, b2::Rect{N}) where N
    for i = 1:N
        maximum(b1)[i] < maximum(b2)[i] && minimum(b1)[i] > minimum(b2)[i] ||
        return false
    end
    true
end

function finishes(b1::Rect{N}, b2::Rect{N}) where N
    if maximum(b1) == maximum(b2)
        for i = 1:N
            minimum(b1)[i] > minimum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

#
# Containment
#
function isinside(rect::SimpleRectangle, x::Real, y::Real)
    rect.x <= x && rect.y <= y && rect.x + rect.w >= x && rect.y + rect.h >= y
end

# TODO only have point in c and deprecate above methods
in(x::AbstractPoint{2}, c::Circle) = isinside(c, x...)
in(x::AbstractPoint{2}, c::SimpleRectangle) = isinside(c, x...)




"""
Check if HyperRectangles are contained in each other. This does not use
strict inequality, so HyperRectangles may share faces and this will still
return true.
"""
function contains(b1::Rect{N}, b2::Rect{N}) where N
    for i = 1:N
        maximum(b2)[i] <= maximum(b1)[i] &&
            minimum(b2)[i] >= minimum(b1)[i] ||
            return false
    end
    return true
end

"""
Check if a point is contained in a Rect. This will return true if
the point is on a face of the Rect.
"""
function contains(b1::Rect{N, T}, pt::Union{FixedVector, AbstractVector}) where {T, N}
    for i = 1:N
        pt[i] <= maximum(b1)[i] && pt[i] >= minimum(b1)[i] || return false
    end
    return true
end

"""
Check if HyperRectangles are contained in each other. This does not use
strict inequality, so HyperRectangles may share faces and this will still
return true.
"""
in(b1::Rect, b2::Rect) = contains(b2, b1)

"""
Check if a point is contained in a Rect. This will return true if
the point is on a face of the Rect.
"""
in(pt::Union{FixedVector, AbstractVector}, b1::Rect) = contains(b1, pt)



#
# Equality
#
(==)(b1::Rect, b2::Rect) = minimum(b1) == minimum(b2) && widths(b1) == widths(b2)


isequal(b1::Rect, b2::Rect) = b1 == b2

isless(a::SimpleRectangle, b::SimpleRectangle) = isless(area(a), area(b))


centered(R::Type{HyperRectangle{N,T}}) where {N, T} = R(Vec{N,T}(-0.5), Vec{N,T}(1))
centered(::Type{T}) where {T <: HyperRectangle} = centered(HyperRectangle{ndims_or(T, 3), eltype_or(T, Float32)})
