"""
    Cone{T}(origin::Point3, tip::Point3, radius)

A Cone is a cylinder where one end has a radius of 0. It is defined by an
`origin` with a finite `radius` which linearly decreases to 0 at the `tip`.
"""
struct Cone{T} <: GeometryPrimitive{3, T}
    origin::Point3{T}
    tip::Point3{T}
    radius::T
end

function Cone(origin::Point3{T1}, tip::Point3{T2}, radius::T3) where {T1, T2, T3}
    T = promote_type(T1, T2, T3)
    return Cone{T}(origin, tip, radius)
end

origin(c::Cone) = c.origin
extremity(c::Cone) = c.tip
radius(c::Cone) = c.radius
height(c::Cone) = norm(c.tip - c.origin)
direction(c::Cone) = (c.tip .- c.origin) ./ height(c)

# Note:
# nvertices is matched with Cylinder, where each end has half the vertices. That
# results in less than nvertices for Cone, but allows a Cylinder and a Cone to
# be seamless matched with the same `nvertices`

function coordinates(c::Cone{T}, nvertices=30) where {T}
    nvertices += isodd(nvertices)
    nhalf = div(nvertices, 2)

    R = cylinder_rotation_matrix(direction(c))
    step = 2pi / nhalf

    ps = Vector{Point3{T}}(undef, nhalf + 2)
    for i in 1:nhalf
        phi = (i-1) * step
        ps[i] = R * Point3{T}(c.radius * cos(phi), c.radius * sin(phi), 0) + c.origin
    end
    ps[end-1] = c.tip
    ps[end] = c.origin

    return ps
end

function normals(c::Cone, nvertices = 30)
    nvertices += isodd(nvertices)
    nhalf = div(nvertices, 2)

    R = cylinder_rotation_matrix(direction(c))
    step = 2pi / nhalf

    ns = Vector{Vec3f}(undef, nhalf + 2)
    # shell at origin
    # normals are angled in z direction due to change in radius (from radius to 0)
    # This can be calculated from triangles
    z = radius(c) / height(c)
    norm = 1.0 / sqrt(1 + z*z)
    for i in 1:nhalf
        phi = (i-1) * step
        ns[i] = R * (norm * Vec3f(cos(phi), sin(phi), z))
    end

    # tip - this is undefined / should be all ring angles at once
    # for rendering it is useful to define this as Vec3f(0), because tip normal
    # has no useful value to contribute to the interpolated fragment normal
    ns[end-1] = Vec3f(0)

    # cap
    ns[end] = Vec3f(normalize(c.origin - c.tip))

    faces = Vector{GLTriangleFace}(undef, nvertices)

    # shell
    for i in 1:nhalf
        faces[i] = GLTriangleFace(i, mod1(i+1, nhalf), nhalf+1)
    end

    # cap
    for i in 1:nhalf
        faces[i+nhalf] = GLTriangleFace(nhalf + 2)
    end

    return FaceView(ns, faces)
end

function faces(::Cone, facets=30)
    nvertices = facets + isodd(facets)
    nhalf = div(nvertices, 2)

    faces = Vector{GLTriangleFace}(undef, nvertices)

    # shell
    for i in 1:nhalf
        faces[i] = GLTriangleFace(i, mod1(i+1, nhalf), nhalf+1)
    end

    # cap
    for i in 1:nhalf
        faces[i+nhalf] = GLTriangleFace(i, mod1(i+1, nhalf), nhalf+2)
    end

    return faces
end
