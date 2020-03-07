function Rect(geometry::AbstractArray{<: Point{N, T}}) where {N,T}
    return Rect{N,T}(geometry)
end

"""
Construct a HyperRectangle enclosing all points.
"""
function (t::Type{Rect{N1, T1}})(
        geometry::AbstractArray{PT}
    ) where {N1, T1, PT <: Point}
    N2, T2 = length(PT), eltype(PT)
    @assert N1 >= N2
    vmin = Point{N2, T2}(typemax(T2))
    vmax = Point{N2, T2}(typemin(T2))
    for p in geometry
        vmin, vmax = minmax(p, vmin, vmax)
    end
    o = vmin
    w = vmax - vmin
    if N1 > N2
        z = zero(Vec{N1-N2, T1})
        return Rect{N1, T1}(vcat(o, z),
                                     vcat(w, z))
    else
        return Rect{N1, T1}(o, w)
   end
end


function Rect{T}(a::Pyramid) where T
    w,h = a.width/T(2), a.length
    m = Vec{3,T}(a.middle)
    Rect{T}(m-Vec{3,T}(w,w,0), m+Vec{3,T}(w, w, h))
end

Rect{T}(a::Cube) where T = Rect{T}(origin(a), widths(a))

Rect{T}(a::AbstractMesh) where T = Rect{T}(vertices(a))
Rect{N, T}(a::AbstractMesh) where {N, T} = Rect{N, T}(vertices(a))
