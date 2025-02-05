# Boundingbox-like Rect constructors

Rect(p::AbstractGeometry{N, T}) where {N, T} = Rect{N, T}(p)
RectT{T}(p::AbstractGeometry{N}) where {N, T} = Rect{N, T}(p)
Rect{N}(p::AbstractGeometry{_N, T}) where {N <: Val, _N, T} = Rect{N, T}(p)

Rect(p::AbstractArray{<: VecTypes{N, T}}) where {N, T} = Rect{N, T}(p)
RectT{T}(p::AbstractArray{<: VecTypes{N}}) where {N, T} = Rect{N, T}(p)
Rect{N}(p::AbstractArray{<: VecTypes{_N, T}}) where {N, _N, T} = Rect{N, T}(p)

# Implementations
# Specialize fully typed Rect constructors
function Rect{N, T}(geom::AbstractGeometry) where {N, T <: Number}
    if applicable(Rect{T}, geom)
        @warn "`Rect{T}(geom)` is deprecated as the final boundingbox method. Define `Rect{N, T}(geom)` instead."
        return Rect{T}(geom)
    else
        return Rect{N, T}(coordinates(geom))
    end
end

function bbox_dim_check(trg, src::Integer)
    @assert trg isa Integer "Rect{$trg, $T1} is invalid. This may have happened due to calling Rect{$N1}(obj) to get a bounding box."
    if trg < src
        throw(ArgumentError("Cannot construct a $trg dimensional bounding box from $src dimensional Points. ($trg must be ≥ $src)"))
    end
end

"""
    Rect(points::AbstractArray{<: VecTypes})

Construct a bounding box containing all the given points.
"""
function Rect{N1, T1}(points::AbstractArray{<: VecTypes{N2, T2}}) where {N1, T1, N2, T2}
    bbox_dim_check(N1, N2)
    vmin = Point{N2,T2}(typemax(T2))
    vmax = Point{N2,T2}(typemin(T2))
    for p in points
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
function Rect{N, T}(a::Pyramid) where {N, T}
    bbox_dim_check(N, 3)
    w, h = a.width, a.length
    m = Vec{3,T}(a.middle)
    return Rect{N, T}(m .- Vec{3,T}(w / T(2), w / T(2), 0), Vec{3,T}(w, w, h))
end

function Rect{N, T}(a::HyperSphere{N2}) where {N, N2, T}
    bbox_dim_check(N, N2)
    return Rect{N, T}(minimum(a), widths(a))
end

# TODO: exact implementation that doesn't rely on coordinates
# function Rect{N, T}(a::Cylinder) where {N, T}
#     return Rect{N, T}(...)
# end