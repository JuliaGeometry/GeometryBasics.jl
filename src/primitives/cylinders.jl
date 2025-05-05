"""
    Cylinder{T}(origin::Point3, extremity::Point3, radius)

A `Cylinder` is a 3D primitive defined by an `origin`, an `extremity` (end point)
and a `radius`.
"""
struct Cylinder{T} <: GeometryPrimitive{3, T}
    origin::Point3{T}
    extremity::Point3{T}
    r::T
end

function Cylinder(origin::Point3{T1}, extremity::Point3{T2}, radius::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    Cylinder{T}(origin, extremity, radius)
end

origin(c::Cylinder) = c.origin
extremity(c::Cylinder) = c.extremity
radius(c::Cylinder) = c.r
height(c::Cylinder) = norm(c.extremity - c.origin)
direction(c::Cylinder) = (c.extremity .- c.origin) ./ height(c)

"""
    cylinder_rotation_matrix(direction::VecTypes{3})

Creates a basis transformation matrix `R` that maps the third dimension to the
given `direction` and the first and second to orthogonal directions. This allows
you to encode a rotation around `direction` in the first two components and
transform it with `R * rotated_point`.
"""
function cylinder_rotation_matrix(d3::VecTypes{3, T}) where {T}
    u = Vec{3, T}(d3[1], d3[2], d3[3])
    if abs(u[1]) > 0 || abs(u[2]) > 0
        v = Vec{3, T}(u[2], -u[1], T(0))
    else
        v = Vec{3, T}(T(0), -u[3], u[2])
    end
    v = normalize(v)
    w = Vec{3, T}(u[2] * v[3] - u[3] * v[2], -u[1] * v[3] + u[3] * v[1],
                  u[1] * v[2] - u[2] * v[1])
    return Mat{3, 3, T}(v..., w..., u...)
end

function coordinates(c::Cylinder{T}, nvertices=30) where {T}
    nvertices += isodd(nvertices)
    nhalf = div(nvertices, 2)

    R = cylinder_rotation_matrix(direction(c))
    step = 2pi / nhalf

    ps = Vector{Point3{T}}(undef, nvertices + 2)
    for i in 1:nhalf
        phi = (i-1) * step
        ps[i] = R * Point3{T}(c.r * cos(phi), c.r * sin(phi), 0) + c.origin
    end
    for i in 1:nhalf
        phi = (i-1) * step
        ps[i + nhalf] = R * Point3{T}(c.r * cos(phi), c.r * sin(phi), 0) + c.extremity
    end
    ps[end-1] = c.origin
    ps[end] = c.extremity

    return ps
end

function normals(c::Cylinder, nvertices = 30)
    nvertices += isodd(nvertices)
    nhalf = div(nvertices, 2)

    R = cylinder_rotation_matrix(direction(c))
    step = 2pi / nhalf

    ns = Vector{Vec3f}(undef, nhalf + 2)
    for i in 1:nhalf
        phi = (i-1) * step
        ns[i] = R * Vec3f(cos(phi), sin(phi), 0)
    end
    ns[end-1] = R * Vec3f(0, 0, -1)
    ns[end] = R * Vec3f(0, 0, 1)

    disk1 = map(i -> GLTriangleFace(nhalf+1), 1:nhalf)
    mantle = map(i -> QuadFace(i, mod1(i+1, nhalf), mod1(i+1, nhalf), i), 1:nhalf)
    disk2 = map(i -> GLTriangleFace(nhalf+2), 1:nhalf)
    fs = vcat(disk1, mantle, disk2)

    return FaceView(ns, fs)
end

function faces(::Cylinder, facets=30)
    nvertices = facets + isodd(facets)
    nhalf = div(nvertices, 2)

    disk1 = map(i -> GLTriangleFace(nvertices+1, mod1(i+1, nhalf), i), 1:nhalf)
    mantle = map(1:nhalf) do i
        i1 = mod1(i+1, nhalf)
        QuadFace(i, i1, i1 + nhalf, i+nhalf)
    end
    disk2 = map(i -> GLTriangleFace(nvertices+2, i+nhalf, mod1(i+1, nhalf)+nhalf), 1:nhalf)

    return vcat(disk1, mantle, disk2)
end
