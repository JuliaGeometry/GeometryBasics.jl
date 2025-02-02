"""
    Pyramid(middle::Point3, length::Real, width::Real)

A Pyramid is an axis-aligned primitive where the tip of the Pyramid extends by 
`length` from `middle` in z direction and the square base extends by `width`
in ±x and ±y direction from `middle`. 
"""
struct Pyramid{T} <: GeometryPrimitive{3,T}
    middle::Point{3,T}
    length::T
    width::T
end

function Pyramid(middle::Point{3, T1}, length::T2, width::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    return Pyramid(Point3{T}(middle), T(length), T(width))
end

function coordinates(p::Pyramid{T}) where {T}
    leftup = Point{3,T}(-p.width, p.width, 0) / 2
    leftdown = Point(-p.width, -p.width, 0) / 2
    tip = Point{3,T}(p.middle + Point{3,T}(0, 0, p.length))
    lu = Point{3,T}(p.middle + leftup)
    ld = Point{3,T}(p.middle + leftdown)
    ru = Point{3,T}(p.middle - leftdown)
    rd = Point{3,T}(p.middle - leftup)
    return Point{3,T}[tip, rd, ru, lu, ld]
end

function normals(p::Pyramid)
    w = p.width; h = p.length
    ns = normalize.(Vec3f[(h, 0, w), (0, h, w), (-h, 0, w), (0, -h, w), (0, 0, -1)])
    fs = [GLTriangleFace(1), GLTriangleFace(2), GLTriangleFace(3), GLTriangleFace(4), QuadFace(5)]
    return FaceView(ns, fs)
end

function faces(::Pyramid)
    return [
        GLTriangleFace(1, 2, 3),
        GLTriangleFace(1, 3, 4),
        GLTriangleFace(1, 4, 5),
        GLTriangleFace(1, 5, 2),
        QuadFace(2, 3, 4, 5)
    ]
end
