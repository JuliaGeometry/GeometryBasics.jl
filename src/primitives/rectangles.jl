
"""
    HyperRectangle{N, T}

A `HyperRectangle` is a generalization of a rectangle into N-dimensions.
Formally it is the Cartesian product of intervals, which is represented by the
`origin` and `widths` fields, whose indices correspond to each of the `N` axes.
"""
struct HyperRectangle{N,T} <: GeometryPrimitive{N,T}
    origin::Vec{N,T}
    widths::Vec{N,T}

    function HyperRectangle{N, T}(o::VecTypes, w::VecTypes) where {N, T}
        return new{N, T}(convert(Vec{N, T}, o), convert(Vec{N, T}, w))
    end
    function HyperRectangle{N, T}(o::VecTypes, w::VecTypes) where {N, T <: Integer}
        return new{N, T}(convert(Vec{N, T}, round.(T, o)), convert(Vec{N, T}, round.(T, w)))
    end
end

##
# Constructors & typealiases

"""
    const Rect{N,T} = HyperRectangle{N,T}

A rectangle in N dimensions, formally the Cartesian product of intervals. See also [`HyperRectangle`](@ref). Its aliases are

|        |`T`(eltype)|`Float64` |`Float32` |`Int`     |
|--------|-----------|----------|----------|----------|
|`N`(dim)|`Rect{N,T}`|`Rectd{N}`|`Rectf{N}`|`Recti{N}`|
|`2`     |`Rect2{T}` |`Rect2d`  |`Rect2f`  |`Rect2i`  |
|`3`     |`Rect3{T}` |`Rect3d`  |`Rect3f`  |`Rect3i`  |

There is an additional unexported alias `RectT` that simply reverses the order of type parameters: `RectT{T,N} == Rect{N,T}`.

"""
Rect, Rect2, Rect3, RectT, Rectd, Rect2d, Rect3d, Rectf, Rect2f, Rect3f, Recti, Rect2i, Rect3i

const Rect{N,T} = HyperRectangle{N,T}
const Rect2{T} = Rect{2,T}
const Rect3{T} = Rect{3,T}

const RectT{T,N} = Rect{N,T}

const Rectd{N} = Rect{N,Float64}
const Rect2d = Rect2{Float64}
const Rect3d = Rect3{Float64}

const Rectf{N} = Rect{N,Float32}
const Rect2f = Rect2{Float32}
const Rect3f = Rect3{Float32}

const Recti{N} = Rect{N,Int}
const Rect2i = Rect2{Int}
const Rect3i = Rect3{Int}


# Constructors


Rect() = Rect{2,Float32}()
RectT{T}() where {T} = Rect{2,T}()
Rect{N}() where {N} = Rect{N,Float32}()

function Rect{N,T}() where {T,N}
    # empty constructor such that update will always include the first point
    return Rect{N,T}(Vec{N,T}(typemax(T)), Vec{N,T}(typemin(T)))
end

# Rect(numbers...)
Rect(args::Vararg{Number, N}) where {N} = Rect{div(N, 2), promote_type(typeof.(args)...)}(args...)
RectT{T}(args::Vararg{Number, N}) where {N, T} = Rect{div(N, 2), T}(args...)
Rect{N}(args::Vararg{Number}) where {N} = Rect{N, promote_type(typeof.(args)...)}(args...)

"""
    Rect(vals::Number...)

Rect constructor for individually specified intervals.
e.g. Rect(0,0,1,2) has origin == Vec(0,0) and
width == Vec(1,2)
"""
function Rect{N, T}(vals::Vararg{Number, M}) where {N, M, T}
    n, r = fldmod(M, 2)
    if r != 0 || n != N
        throw(ArgumentError("Number of arguments must be compatible with given or derived Rect size. Got $M arguments for a Rect{$N} expecting $(2 * N)."))
    end
    origin, widths = ntuple(i -> vals[i], N), ntuple(i -> vals[i+N], N)
    return Rect{N, T}(Vec(origin), Vec(widths))
end

# VecTypes

     Rect(o::VecTypes{N, T1}, w::VecTypes{N, T2}) where {N, T1, T2} = Rect{N, promote_type(T1, T2)}(o, w)
RectT{ T}(o::VecTypes{N},     w::VecTypes{N})     where {N, T}      = Rect{N, T}(o, w)
Rect{N  }(o::VecTypes{N, T1}, w::VecTypes{N, T2}) where {N, T1, T2} = Rect{N, promote_type(T1, T2)}(o, w)

# Arrays

      Rect(o::AbstractVector{T1}, w::AbstractVector{T2}) where {T1, T2}    = RectT{promote_type(T1, T2)}(Vec(o...), Vec(w...))
RectT{  T}(o::AbstractVector,     w::AbstractVector)     where {T}         = RectT{T}(Vec(o...), Vec(w...))
Rect{N   }(o::AbstractVector{T1}, w::AbstractVector{T2}) where {N, T1, T2} = Rect{N, promote_type(T1, T2)}(Vec(o...), Vec(w...))
Rect{N, T}(o::AbstractVector, w::AbstractVector) where {N, T}              = Rect{N, T}(Vec(o...), Vec(w...))

# mixed number - vectype

      Rect(o::VecTypes{N, <:Number}, args::Vararg{Number, N}) where {N}    = Rect{N   }(o, promote(args...))
RectT{  T}(o::VecTypes{N, <:Number}, args::Vararg{Number, N}) where {N, T} = Rect{N, T}(o, promote(args...))
Rect{N   }(o::VecTypes{N, <:Number}, args::Vararg{Number, N}) where {N}    = Rect{N   }(o, promote(args...))
Rect{N, T}(o::VecTypes{N, <:Number}, args::Vararg{Number, N}) where {N, T} = Rect{N, T}(o, promote(args...))

      Rect(x::Number, y::Number, w::VecTypes{2, <:Number})           = Rect{2   }(Vec(x, y), w)
RectT{  T}(x::Number, y::Number, w::VecTypes{2, <:Number}) where {T} = Rect{2, T}(Vec(x, y), w)
Rect{2   }(x::Number, y::Number, w::VecTypes{2, <:Number})           = Rect{2   }(Vec(x, y), w)
Rect{2, T}(x::Number, y::Number, w::VecTypes{2, <:Number}) where {T} = Rect{2, T}(Vec(x, y), w)

      Rect(x::Number, y::Number, z::Number, w::VecTypes{3, <:Number})           = Rect{3   }(Vec(x, y, z), w)
RectT{  T}(x::Number, y::Number, z::Number, w::VecTypes{3, <:Number}) where {T} = Rect{3, T}(Vec(x, y, z), w)
Rect{3   }(x::Number, y::Number, z::Number, w::VecTypes{3, <:Number})           = Rect{3   }(Vec(x, y, z), w)
Rect{3, T}(x::Number, y::Number, z::Number, w::VecTypes{3, <:Number}) where {T} = Rect{3, T}(Vec(x, y, z), w)

# copy constructors (allow explicit truncation)

      Rect(r::Rect{N, T}) where {N, T}      = Rect{N, T}(origin(r), widths(r))
RectT{  T}(r::Rect{N})    where {N, T}      = Rect{N, T}(origin(r), widths(r))
Rect{N   }(r::Rect{N2, T}) where {N, N2, T} = Rect{N, T}(Vec{min(N, N2), T}(origin(r)), Vec{min(N, N2), T}(widths(r)))
Rect{N, T}(r::Rect{N2})    where {N, N2, T} = Rect{N, T}(Vec{min(N, N2), T}(origin(r)), Vec{min(N, N2), T}(widths(r)))

# dimensional promotion

Rect{3, T}(o::VecTypes{2}, w::VecTypes{3}) where {T} = Rect{3, T}(Vec(o..., 0), w)
Rect{3, T}(o::VecTypes{3}, w::VecTypes{2}) where {T} = Rect{3, T}(o,            Vec(w..., 0))
Rect{3, T}(o::VecTypes{2}, w::VecTypes{2}) where {T} = Rect{3, T}(Vec(o..., 0), Vec(w..., 0))

# centered Rects

centered(R::Type{Rect{N,T}}) where {N,T} = R(Vec{N,T}(-0.5), Vec{N,T}(1))
centered(R::Type{RectT{T}}) where {T} = R(Vec{2,T}(-0.5), Vec{2,T}(1))
centered(R::Type{Rect{N}}) where {N} = R(Vec{N,Float32}(-0.5), Vec{N,Float32}(1))
centered(R::Type{Rect}) = R(Vec{2,Float32}(-0.5), Vec{2,Float32}(1))

# TODO These are kinda silly
function Rect2(xy::NamedTuple{(:x, :y)}, wh::NamedTuple{(:width, :height)})
    return Rect2(xy.x, xy.y, wh.width, wh.height)
end

Rect(ow::Tuple) = Rect(ow...)
RectT{T}(ow::Tuple) where {T} = RectT{T}(ow...)
Rect{N}(ow::Tuple) where {N} = Rect{N}(ow...)
Rect{N, T}(ow::Tuple) where {N, T} = Rect{N, T}(ow...)


# Utilities


# allow auto-conversion between different eltypes
Base.convert(::Type{Rect{N, T}}, r::Rect{N}) where {N, T} = Rect{N, T}(r)

origin(prim::Rect) = prim.origin
Base.maximum(prim::Rect) = origin(prim) + widths(prim)
Base.minimum(prim::Rect) = origin(prim)
Base.length(prim::Rect{N,T}) where {T,N} = N
widths(prim::Rect) = prim.widths

width(prim::Rect) = prim.widths[1]
height(prim::Rect) = prim.widths[2]

volume(prim::HyperRectangle) = prod(prim.widths)
area(prim::Rect2) = volume(prim)

# function Base.round(::Type{Rect{N, T}}, x::Rect{N}, mode::RoundingMode=RoundNearest) where {N, T}
#     mini = round.(T, minimum(x))
#     maxi = round.(T, maximum(x))
#     return Rect{N, T}(mini, maxi .- mini)
# end

"""
    split(rectangle, axis, value)

Splits an Rect into two along an axis at a given location.
"""
split(b::Rect, axis, value::Integer) = _split(b, axis, value)
split(b::Rect, axis, value::Number) = _split(b, axis, value)
function _split(b::H, axis, value) where {H<:Rect}
    bmin = minimum(b)
    bmax = maximum(b)
    b1max = Base.setindex(bmax, value, axis)
    b2min = Base.setindex(bmin, value, axis)

    return H(bmin, b1max - bmin), H(b2min, bmax - b2min)
end

###
# Algebraic operations

"""
    *(m::Mat, h::Rect)

Transform a `Rect` using a matrix. Maintains axis-align properties
so a significantly larger Rect may be generated.
"""
function Base.:(*)(m::Mat{N1,N1,T1}, h::Rect{N2,T2}) where {N1,N2,T1,T2}

    # TypeVar constants
    T = promote_type(T1, T2)
    D = N1 - N2

    # get all points on the Rect
    d = decompose(Point, h)
    # make sure our points are sized for the transform
    pts = (Vec{N1,T}[vcat(pt, ones(Vec{D,T})) for pt in d]...,)::NTuple{2^N2,Vec{N1,T}}

    vmin = Vec{N1,T}(typemax(T))
    vmax = Vec{N1,T}(typemin(T))
    # transform all points, tracking min and max points
    for pt in pts
        pn = m * pt
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    return Rect{N2,T}(vmin, vmax - vmin)
end

# equal dimension case
function Base.:(*)(m::Mat{N,N,T1}, h::Rect{N,T2}) where {N,T1,T2}

    # TypeVar constants
    T = promote_type(T1, T2)

    # get all points on the Rect
    pts = decompose(Point, h)

    # make sure our points are sized for the transform
    vmin = Vec{N,T}(typemax(T))
    vmax = Vec{N,T}(typemin(T))

    # transform all points, tracking min and max points
    for pt in pts
        pn = m * Vec(pt)
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    return Rect{N,T}(vmin, vmax - vmin)
end

# fast path. TODO make other versions fast without code duplications like now
function Base.:(*)(m::Mat{4,4,T}, h::Rect{3,T}) where {T}
    # equal dimension case

    # get all points on the Rect
    pts = (Vec{4,T}(0.0, 0.0, 0.0, 1.0), Vec{4,T}(1.0, 0.0, 0.0, 1.0),
           Vec{4,T}(0.0, 1.0, 0.0, 1.0), Vec{4,T}(1.0, 1.0, 0.0, 1.0),
           Vec{4,T}(0.0, 0.0, 1.0, 1.0), Vec{4,T}(1.0, 0.0, 1.0, 1.0),
           Vec{4,T}(0.0, 1.0, 1.0, 1.0), Vec{4,T}(1.0, 1.0, 1.0, 1.0))

    # make sure our points are sized for the transform
    vmin = Vec{4,T}(typemax(T))
    vmax = Vec{4,T}(typemin(T))
    o, w = origin(h), widths(h)
    _o = Vec{4,T}(o[1], o[2], o[3], T(0))
    _w = Vec{4,T}(w[1], w[2], w[3], T(1))
    # transform all points, tracking min and max points
    for pt in pts
        pn = m * (_o + (pt .* _w))
        vmin = min.(pn, vmin)
        vmax = max.(pn, vmax)
    end
    _vmin = Vec{3,T}(vmin[1], vmin[2], vmin[3])
    _vmax = Vec{3,T}(vmax[1], vmax[2], vmax[3])
    return Rect{3,T}(_vmin, _vmax - _vmin)
end

Base.:(-)(h::Rect{N,T}, move::Number) where {N,T} = h - Vec{N,T}(move)
Base.:(+)(h::Rect{N,T}, move::Number) where {N,T} = h + Vec{N,T}(move)

function Base.:(-)(h::Rect{N,T}, move::StaticVector{N}) where {N,T}
    return Rect{N,T}(minimum(h) .- move, widths(h))
end

function Base.:(+)(h::Rect{N,T}, move::StaticVector{N}) where {N,T}
    return Rect{N,T}(minimum(h) .+ move, widths(h))
end

function Base.:(*)(rect::Rect, scaling::Union{Number,StaticVector})
    return Rect(minimum(rect) .* scaling, widths(rect) .* scaling)
end

# Enables rectangular indexing into a matrix
function Base.to_indices(A::AbstractMatrix{T}, I::Tuple{Rect2{IT}}) where {T,IT<:Integer}
    rect = I[1]
    mini = minimum(rect)
    wh = widths(rect)
    return ((mini[1] + 1):(mini[1] + wh[1]), (mini[2] + 1):(mini[2] + wh[2]))
end

function _minmax(p::StaticVector, vmin, vmax)
    any(isnan, p) && return (vmin, vmax)
    return min.(p, vmin), max.(p, vmax)
end

# TODO: doesn't work regardless
# Annoying special case for view(Vector{Point}, Vector{Face})
# function Base.minmax(tup::Tuple, vmin, vmax)
#     for p in tup
#         any(isnan, p) && continue
#         vmin = min.(p, vmin)
#         vmax = max.(p, vmax)
#     end
#     return vmin, vmax
# end

function positive_widths(rect::Rect{N,T}) where {N,T}
    mini, maxi = minimum(rect), maximum(rect)
    realmin = min.(mini, maxi)
    realmax = max.(mini, maxi)
    return Rect{N,T}(realmin, realmax .- realmin)
end

###
# set operations

"""
    isempty(h::Rect)

Return `true` if any of the widths of `h` are negative.
"""
Base.isempty(h::Rect{N,T}) where {N,T} = any(<(zero(T)), h.widths)

"""
    union(r1::Rect{N}, r2::Rect{N})

Returns a new `Rect{N}` which contains both r1 and r2.
"""
function Base.union(h1::Rect{N}, h2::Rect{N}) where {N}
    m = min.(minimum(h1), minimum(h2))
    mm = max.(maximum(h1), maximum(h2))
    return Rect{N}(m, mm - m)
end

# TODO: What should this be? The difference is "h2 - h1", which could leave an
# L shaped cutout. Should we pad that back out into a full rect?
# """
#     diff(h1::Rect, h2::Rect)

# Perform a difference between two Rects.
# """
# diff(h1::Rect, h2::Rect) = h1

"""
    intersect(h1::Rect, h2::Rect)

Perform a intersection between two Rects.
"""
function Base.intersect(h1::Rect{N}, h2::Rect{N}) where {N}
    m = max.(minimum(h1), minimum(h2))
    mm = min.(maximum(h1), maximum(h2))
    return Rect{N}(m, mm - m)
end

function update(b::Rect{N,T}, v::VecTypes{N,T2}) where {N,T,T2}
    return update(b, Vec{N,T}(v))
end

function update(b::Rect{N,T}, v::VecTypes{N,T}) where {N,T}
    m = min.(minimum(b), v)
    maxi = maximum(b)
    mm = if any(isnan, maxi)
        v - m
    else
        max.(v, maxi) - m
    end
    return Rect{N,T}(m, mm)
end

# Min maximum distance functions between hrectangle and point for a given dimension
function min_dist_dim(rect::Rect{N,T}, p::VecTypes{N,T}, dim::Int) where {N,T}
    return max(zero(T), max(minimum(rect)[dim] - p[dim], p[dim] - maximum(rect)[dim]))
end

function max_dist_dim(rect::Rect{N,T}, p::VecTypes{N,T}, dim::Int) where {N,T}
    return max(maximum(rect)[dim] - p[dim], p[dim] - minimum(rect)[dim])
end

function min_dist_dim(rect1::Rect{N,T}, rect2::Rect{N,T}, dim::Int) where {N,T}
    return max(zero(T),
               max(minimum(rect1)[dim] - maximum(rect2)[dim],
                   minimum(rect2)[dim] - maximum(rect1)[dim]))
end

function max_dist_dim(rect1::Rect{N,T}, rect2::Rect{N,T}, dim::Int) where {N,T}
    return max(maximum(rect1)[dim] - minimum(rect2)[dim],
               maximum(rect2)[dim] - minimum(rect1)[dim])
end

# Total minimum maximum distance functions
function min_euclideansq(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    minimum_dist = T(0.0)
    for dim in 1:length(p)
        d = min_dist_dim(rect, p, dim)
        minimum_dist += d * d
    end
    return minimum_dist
end

function max_euclideansq(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    maximum_dist = T(0.0)
    for dim in 1:length(p)
        d = max_dist_dim(rect, p, dim)
        maximum_dist += d * d
    end
    return maximum_dist
end

function min_euclidean(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    return sqrt(min_euclideansq(rect, p))
end

function max_euclidean(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    return sqrt(max_euclideansq(rect, p))
end

# Functions that return both minimum and maximum for convenience
function minmax_dist_dim(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}},
                         dim::Int) where {N,T}
    minimum_d = min_dist_dim(rect, p, dim)
    maximum_d = max_dist_dim(rect, p, dim)
    return minimum_d, maximum_d
end

function minmax_euclideansq(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    minimum_dist = min_euclideansq(rect, p)
    maximum_dist = max_euclideansq(rect, p)
    return minimum_dist, maximum_dist
end

function minmax_euclidean(rect::Rect{N,T}, p::Union{VecTypes{N,T},Rect{N,T}}) where {N,T}
    minimumsq, maximumsq = minmax_euclideansq(rect, p)
    return sqrt(minimumsq), sqrt(maximumsq)
end

# http://en.wikipedia.org/wiki/Allen%27s_interval_algebra
function before(b1::Rect{N}, b2::Rect{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

meets(b1::Rect{N}, b2::Rect{N}) where {N} = maximum(b1) == minimum(b2)

function overlaps(b1::Rect{N}, b2::Rect{N}) where {N}
    for i in 1:N
        maximum(b2)[i] > maximum(b1)[i] > minimum(b2)[i] &&
            minimum(b1)[i] < minimum(b2)[i] || return false
    end
    return true
end

function starts(b1::Rect{N}, b2::Rect{N}) where {N}
    return if minimum(b1) == minimum(b2)
        for i in 1:N
            maximum(b1)[i] < maximum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

function during(b1::Rect{N}, b2::Rect{N}) where {N}
    for i in 1:N
        maximum(b1)[i] < maximum(b2)[i] && minimum(b1)[i] > minimum(b2)[i] || return false
    end
    return true
end

function finishes(b1::Rect{N}, b2::Rect{N}) where {N}
    return if maximum(b1) == maximum(b2)
        for i in 1:N
            minimum(b1)[i] > minimum(b2)[i] || return false
        end
        return true
    else
        return false
    end
end

###
# Containment

"""
    in(b1::Rect, b2::Rect)

Check if Rect `b1` is contained in `b2`. This does not use
strict inequality, so Rects may share faces and this will still
return true.
"""
function Base.in(b1::Rect{N}, b2::Rect{N}) where {N}
    for i in 1:N
        maximum(b1)[i] <= maximum(b2)[i] && minimum(b1)[i] >= minimum(b2)[i] || return false
    end
    return true
end

"""
    in(pt::VecTypes, b1::Rect{N, T})

Check if a point is contained in a Rect. This will return true if
the point is on a face of the Rect.
"""
function Base.in(pt::VecTypes, b1::Rect{N,T}) where {T,N}
    for i in 1:N
        pt[i] <= maximum(b1)[i] && pt[i] >= minimum(b1)[i] || return false
    end
    return true
end

#
# Equality
#
Base.:(==)(b1::Rect, b2::Rect) = origin(b1) == origin(b2) && widths(b1) == widths(b2)

function Base.isapprox(r1::Rect, r2::Rect; kwargs...)
    return isapprox(origin(r1), origin(r2); kwargs...) && isapprox(widths(r1), widths(r2); kwargs...)
end

##
# Rect2 decomposition

function faces(rect::Rect2, nvertices=(2, 2))
    if nvertices == (2, 2)
        return [QuadFace(1,2,3,4)]
    else
        w, h = nvertices
        idx = LinearIndices(nvertices)
        quad(i, j) = QuadFace{Int}(idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1])
        return [quad(i, j) for j in 1:(h - 1) for i in 1:(w - 1)]
    end
end

function coordinates(rect::Rect2{T}, nvertices=(2, 2)) where {T}
    mini, maxi = extrema(rect)
    if nvertices == (2, 2)
        return Point2{T}[mini, (maxi[1], mini[2]), maxi, (mini[1], maxi[2])]
    else
        xrange, yrange = LinRange.(mini, maxi, nvertices)
        return [Point(x, y) for y in yrange for x in xrange]
    end
end

function texturecoordinates(rect::Rect2{T}, nvertices=(2, 2)) where {T}
    ps = coordinates(Rect2{T}(0, 0, 1, 1), nvertices)
    ps = [Vec2{T}(0, 1) .+ Vec2{T}(1, -1) .* p for p in ps]
    return ps
end

function normals(rect::Rect2, nvertices=(2, 2))
    return Iterators.repeated((0, 0, 1), prod(nvertices))
end

##
# Rect3 decomposition
function coordinates(rect::Rect3{T}) where T
    # TODO use n
    w = widths(rect)
    o = origin(rect)
    return Point{3, T}[o + (x, y, z) .* w for x in (0, 1) for y in (0, 1) for z in (0, 1)]
end

function normals(::Rect3)
    ns = Vec3f[(-1,0,0), (1,0,0), (0,-1,0), (0,1,0), (0,0,-1), (0,0,1)]
    return FaceView(ns, QuadFace{Int}.(1:6))
end

function texturecoordinates(rect::Rect3)
    return coordinates(Rect3(0, 0, 0, 1, 1, 1))
end

function faces(::Rect3)
    return QuadFace{Int}[
        (1, 2, 4, 3), (7, 8, 6, 5), (5, 6, 2, 1),
        (3, 4, 8, 7), (1, 3, 7, 5), (6, 8, 4, 2)
    ]
end

function Extents.extent(rect::Rect2)
    (xmin, ymin), (xmax, ymax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end
function Extents.extent(rect::Rect3)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax), Z=(zmin, zmax))
end
