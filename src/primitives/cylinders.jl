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

function rotation(c::Cylinder{T}) where {T}
    d3 = direction(c)
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

    R = rotation(c)
    h = height(c)
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

    R = rotation(c)
    step = 2pi / nhalf
    
    ns = Vector{Vec3f}(undef, nhalf + 2)
    for i in 1:nhalf
        phi = (i-1) * step
        ns[i] = R * Vec3f(cos(phi), sin(phi), 0)
    end
    ns[end-1] = R * Vec3f(0, 0, -1)
    ns[end] = R * Vec3f(0, 0, 1)
    
    return ns
end

function faces(c::Cylinder, facets=30)
    nvertices = facets + isodd(facets)
    nhalf = div(nvertices, 2)

    disk1 = map(1:nhalf) do i
        NormalFace(
            GLTriangleFace(nvertices+1, i, mod1(i+1, nhalf)),
            GLTriangleFace(nhalf+1)
        )
    end

    mantle1 = map(1:nhalf) do i
        i1 = mod1(i+1, nhalf)
        NormalFace(
            GLTriangleFace(i, i+nhalf, i1 + nhalf),
            GLTriangleFace(i, i, i1)
        )
    end

    mantle2 = map(1:nhalf) do i
        i1 = mod1(i+1, nhalf)
        NormalFace(
            GLTriangleFace(i1 + nhalf, i1, i),
            GLTriangleFace(i1, i1, i)
        )
    end

    disk2 = map(1:nhalf) do i
        NormalFace(
            GLTriangleFace(nvertices+2, i+nhalf, mod1(i+1, nhalf)+nhalf),
            GLTriangleFace(nhalf+2)
        )
    end

    return vcat(disk1, mantle1, mantle2, disk2)
end
