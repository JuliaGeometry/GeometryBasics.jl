#=
Functions to compute distances from points to primatives
=#

using LinearAlgebra: normalize,norm,⋅
"""
    closest_point_to_tri(p, a, b, c)

This method from Ericson, "Real-time collision detection"
2005, doesn't require the normal or cross products.
"""
@fastmath function closest_point_to_tri(p,a,b,c)
    # is point `a` closest?
    ab,ac,ap = b-a,c-a,p-a
    d1,d2 = ab ⋅ ap, ac ⋅ ap
    d1 <= 0 && d2 <= 0 && return a

    # is point `b` closest?
    bp = p-b
    d3,d4 = ab ⋅ bp, ac ⋅ bp
    d3 >= 0 && d4 <= d3 && return b

    # is point `c` closest?
    cp = p-c
    d5,d6 = ab ⋅ cp, ac ⋅ cp
    d6 >= 0 && d5 <= d6 && return c

    # is segment `ab` closest?
    vc = d1*d4 - d3*d2
    if (vc <= 0 && d1 >= 0 && d3 <= 0)
        d1==d3 && @show d1,d3
        return a + ab * d1 * inv(d1 - d3)
    end

    # is segment `ac` closest?
    vb = d5*d2 - d1*d6
    if (vb <= 0 && d2 >= 0 && d6 <= 0)
        return a + ac * d2 * inv(d2 - d6)
    end

    # is segment `bc` closest?
    va = d3*d6 - d5*d4
    if (va <= 0 && d4 >= d3 && d5 >= d6)
        return b + (c - b) * (d4 - d3) * inv(d4 - d3 + d5 - d6)
    end

    # closest is interior to `abc`
    denom = inv(va + vb + vc)
    v,w = vb * denom, vc * denom
    return a + ab * v + ac * w
end
"""
    closest(p,tri::Triangle)

Determine the closest point on triangle `tri` to point `p`.
"""
closest(p,tri::Triangle) = closest_point_to_tri(p,tri.points.data...)

"""
    absolute_distance(p,tri::Triangle)

Determine the absolute distance from point `p` to the closest point on triangle `tri`.
"""
absolute_distance(p,tri::Triangle) = norm(p-closest(p,tri))

"""
    signed_distance(p,tri::Triangle)

Determine the signed distance from point `p` to the plane defined by triangle `tri`.
Note that the sign depends on triangle point ordering and `d=0` only implies `p`
lies on the plane, not that it is within the triangle.
"""
signed_distance(p,tri::Triangle) = signed_distance_to_tri(p,tri.points.data...)
signed_distance_to_tri(p,a,b,c) = normalize(orthogonal_vector(a,b,c))⋅(p-a)

nonzero(mesh::AbstractMesh) = (i for i in mesh if sum(abs2,orthogonal_vector(i.points.data...))>0)
"""
    absolute_distance(p, mesh::AbstractMesh)

Return the minimum absolute distance from point `p` to any Triangle in `mesh`.
"""
absolute_distance(p,mesh::AbstractMesh) = minimum(t->absolute_distance(p,t), nonzero(mesh))

"""
    signed_distance(p, mesh::AbstractMesh)

Return the signed distance to the geometry defined by the union of planes
passing through each `mesh` face.
"""
signed_distance(p,mesh::AbstractMesh) = maximum(t->signed_distance(p,t), nonzero(mesh))

"""
    absolute_distance(p,prim) = |signed_distance(p,prim)|

Fallback absolute distance function.
"""
absolute_distance(p,prim) = abs(signed_distance(p,prim))

"""
    signed_distance(p,s::HyperSphere)

Return the signed distance from `p` to `s`.
"""
signed_distance(p,s::HyperSphere) = norm(p-origin(s))-radius(s)

"""
    signed_distance(p,r::Rect)

Return the signed distance from `p` to `r`.
"""
function signed_distance(p,r::Rect)
    # reflect p to positive Rect quadrant and get vector relative to Rect corner
    @show p,origin(r),width(r)
    q = abs.(p-origin(r).-0.5*width(r)).-0.5*width(r)
    # 2-norm for positive distances, ∞-norm for negative
    return norm(max.(q, 0.))+min(maximum(q),0.)
end
