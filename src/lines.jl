
"""
    intersects(a::Line, b::Line) -> Bool, Point

Intersection of 2 line segments `a` and `b`.
Returns `(intersection_found::Bool, intersection_point::Point)`
"""
# 2D Line-segment intersection algorithm by Paul Bourke and many others.
# http://paulbourke.net/geometry/pointlineplane/
function intersects(a::Line{2,T1}, b::Line{2,T2}) where {T1,T2}
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
    if 0 <= unknown_a <= 1 && 0 <= unknown_b <= 1
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
    self_intersections(points::AbstractVector{AbstractPoint})

Finds all self intersections of polygon `points`
"""
function self_intersections(points::AbstractVector{<:AbstractPoint})
    sections = similar(points, 0)
    intersections = Int[]

    wraparound(i) = mod1(i, length(points) - 1)

    for (i, (a, b)) in enumerate(consecutive_pairs(points))
        for (j, (a2, b2)) in enumerate(consecutive_pairs(points))
            is1, is2 = wraparound(i + 1), wraparound(i - 1)
            if i != j &&
               is1 != j &&
               is2 != j &&
               !(i in intersections) &&
               !(j in intersections)
                intersected, p = intersects(Line(a, b), Line(a2, b2))
                if intersected
                    push!(intersections, i, j)
                    push!(sections, p)
                end
            end
        end
    end
    return intersections, sections
end

"""
    split_intersections(points::AbstractVector{AbstractPoint})

Splits polygon `points` into it's self intersecting parts. Only 1 intersection
is handled right now.
"""
function split_intersections(points::AbstractVector{<:AbstractPoint})
    intersections, sections = self_intersections(points)
    return if isempty(intersections)
        return [points]
    elseif length(intersections) == 2 && length(sections) == 1
        a, b = intersections
        p = sections[1]
        a, b = min(a, b), max(a, b)
        poly1 = simple_concat(points, (a + 1):(b - 1), p)
        poly2 = simple_concat(points, (b + 1):(length(points) + a), p)
        return [poly1, poly2]
    else
        error("More than 1 intersections can't be handled currently. Found: $intersections, $sections")
    end
end
