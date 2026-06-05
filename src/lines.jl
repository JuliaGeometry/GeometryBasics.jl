
"""
    intersects(a::Line, b::Line) -> Bool, Point

Intersection of 2 line segments `a` and `b`.
Returns `(intersection_found::Bool, intersection_point::Point)`
"""
# 2D Line-segment intersection algorithm by Paul Bourke and many others.
# http://paulbourke.net/geometry/pointlineplane/
function intersects(a::Line{2,T1}, b::Line{2,T2}; eps = 0) where {T1,T2}
    T = promote_type(T1, T2)
    p0 = zero(Point2{T})

    x1, y1 = a[1]
    x2, y2 = a[2]
    x3, y3 = b[1]
    x4, y4 = b[2]

    denominator = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1))
    numerator_a = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3))
    numerator_b = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3))

    if denominator == 0
        # no intersection: lines are parallel
        return false, p0
    end

    # If we ever need to know if the lines are coincident, we can get that too:
    # denominator == numerator_a == numerator_b == 0 && return :coincident_lines

    # unknown_a and b tell us how far along the line segment the intersection is.
    unknown_a = numerator_a / denominator
    unknown_b = numerator_b / denominator

    # Values between [0, 1] mean the intersection point of the lines rests on
    # both of the line segments.
    if eps <= unknown_a <= 1-eps && eps <= unknown_b <= 1-eps
        # Substituting an unknown back lets us find the intersection point.
        x = x1 + (unknown_a * (x2 - x1))
        y = y1 + (unknown_a * (y2 - y1))
        return true, Point2{T}(x, y)
    end

    # lines intersect, but outside of at least one of these line segments.
    return false, p0
end

function simple_concat(vec::AbstractVector, range, endpoint::P) where {P}
    result = Vector{P}(undef, length(range) + 1)
    for (i, j) in enumerate(range)
        result[i] = vec[mod1(j, length(vec))]
    end
    result[end] = endpoint
    return result
end

function consecutive_pairs(arr)
    n = length(arr)
    return zip(view(arr, 1:(n - 1)), view(arr, 2:n))
end

"""
    self_intersections(points::AbstractVector{<:Point})

Finds all self intersections of in a continuous line described by `points`.
Returns a Vector of indices where each pair `v[2i], v[2i+1]` refers two
intersecting line segments by their first point, and a Vector of intersection
points.

Note that if two points are the same, they will generate a self intersection
unless they are consecutive segments. (The first and last point are assumed to
be shared between the first and last segment.)
"""
function self_intersections(points::AbstractVector{<:VecTypes{D, T}}) where {D, T}
    ti, sections = _self_intersections(points)
    # convert array of tuples to flat array
    return [x for t in ti for x in t], sections
end

function _self_intersections(points::AbstractVector{<:VecTypes{D, T}}) where {D, T}
    sections = similar(points, 0)
    intersections = Tuple{Int, Int}[]

    N = length(points)

    for i in 1:N-3
        a = points[i]; b = points[i+1]
        # i+1 == j describes consecutive segments which are always "intersecting"
        # at point i+1/j. Skip those (start at i+2)
        # Special case: We assume points[1] == points[end] so 1 -> 2 and N-1 -> N
        # always "intersect" at 1/N. Skip this too (end at N-2 in this case)
        for j in i+2 : N-1 - (i == 1)
            a2 = points[j]; b2 = points[j+1]
            intersected, p = intersects(Line(a, b), Line(a2, b2))
            if intersected
                push!(intersections, (i, j))
                push!(sections, p)
            end
        end
    end
    return intersections, sections
end

"""
    split_intersections(points::AbstractVector{<: Point})

Splits polygon `points` into it's self intersecting parts. Only 1 intersection
is handled right now.
"""
function split_intersections(points::AbstractVector{<:VecTypes{N, T}}) where {N, T}
    intersections, sections = _self_intersections(points)
    return if isempty(intersections)
        return [points]
    elseif length(intersections) == 1 && length(sections) == 1
        a, b = intersections[1]
        p = sections[1]
        poly1 = simple_concat(points, (a + 1):b, p)
        poly2 = simple_concat(points, (b + 1):(length(points) + a), p)
        return [poly1, poly2]
    else
        error("More than 1 intersections can't be handled currently. Found: $intersections, $sections")
    end
end
