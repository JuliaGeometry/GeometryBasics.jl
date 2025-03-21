#=
ported from:
http://docs.ros.org/jade/api/convex_decomposition/html/triangulate_8cpp_source.html

The MIT License (MIT)

Copyright (c) 2006 John W. Ratcliff

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
"""
    area(vertices::AbstractVector{Point{3}}, face::NgonFace)

Calculate the area of one face.
"""
function area(vertices::AbstractVector{<:Point}, face::NgonFace)
    return 0.5 * norm(orthogonal_vector(vertices[face]))
end

"""
    area(vertices::AbstractVector{Point}, faces::AbstractVector{NgonFace})

Calculate the area of all faces.
"""
function area(vertices::AbstractVector{<:Point}, faces::AbstractVector{<:NgonFace})
    return sum(x -> area(vertices, x), faces)
end

"""
    area(vertices::AbstractVector{Point}})

Calculate the area of a polygon.

For 2D points, the oriented area is returned (negative when the points are
oriented clockwise).
"""
function area(vertices::AbstractVector{Point{2,T}}) where {T}
    length(vertices) < 3 && return zero(T)
    return T(0.5) * orthogonal_vector(vertices)[3]
end

function area(vertices::AbstractVector{Point{3,T}}) where {T}
    length(vertices) < 3 && return zero(T)
    return T(0.5) * norm(orthogonal_vector(vertices))
end


"""
    in(point, triangle)

Determine if a point is inside of a triangle.
"""
function Base.in(P::T, triangle::Triangle) where {T<:Point}
    A, B, C = coordinates(triangle)

    ap = A .- P
    bp = B .- P
    cp = C .- P

    u = cross(bp, cp)
    v = cross(cp, ap)
    w = cross(ap, bp)

    t0 = zero(eltype(T))
    dot(u, v) < t0 && return false
    dot(u, w) < t0 && return false
    return true
end


function snip(contour::AbstractVector{<:Point}, u, v, w, n, V)
    A = contour[V[u]]
    B = contour[V[v]]
    C = contour[V[w]]
    x = (((B[1] - A[1]) * (C[2] - A[2])) - ((B[2] - A[2]) * (C[1] - A[1])))
    if 0.0000000001f0 > x
        return false
    end

    for p in 1:n
        ((p == u) || (p == v) || (p == w)) && continue
        P = contour[V[p]]
        if P in Triangle(A, B, C)
            return false
        end
    end
    return true
end

"""
    decompose(facetype, contour::AbstractArray{<:Point})

Triangulate a Polygon without hole.

Returns a Vector{`facetype`} defining indexes into `contour`.
"""
function decompose(::Type{F}, points::AbstractVector{<:Point}) where {F<:AbstractFace}
    #= allocate and initialize list of Vertices in polygon =#
    result = F[]

    # the algorithm doesn't like closed contours
    contour = if isapprox(last(points), first(points))
        @view points[1:(end - 1)]
    else
        @view points[1:end]
    end

    n = length(contour)
    if n < 3
        error("Not enough points in the contour. Found: $contour")
    end

    #= we want a counter-clockwise polygon in V =#
    if 0 < area(contour)
        V = Int[i for i in 1:n]
    else
        V = Int[i for i in n:-1:1]
    end

    nv = n

    #=  remove nv-2 Vertices, creating 1 triangle every time =#
    count = 2 * nv   #= error detection =#
    v = nv
    while nv > 2
        #= if we loop, it is probably a non-simple polygon =#
        if 0 >= count
            return result
        end
        count -= 1

        #= three consecutive vertices in current polygon, <u,v,w> =#
        u = v
        (u > nv) && (u = 1) #= previous =#
        v = u + 1
        (v > nv) && (v = 1) #= new v =#
        w = v + 1
        (w > nv) && (w = 1) #= next =#
        if snip(contour, u, v, w, nv, V)
            #= true names of the vertices =#
            a = V[u]
            b = V[v]
            c = V[w]
            #= output Triangle =#
            push!(result, convert_simplex(F, TriangleFace(a, b, c))...)
            #= remove v from remaining polygon =#
            s = v
            t = v + 1
            while t <= nv
                V[s] = V[t]
                s += 1
                t += 1
            end
            nv -= 1
            #= reset error detection counter =#
            count = 2 * nv
        end
    end

    return result
end

function earcut_triangulate(polygon::Vector{Vector{Point{2,Float64}}})
    lengths = map(x -> UInt32(length(x)), polygon)
    len = UInt32(length(lengths))
    array = ccall((:u32_triangulate_f64, libearcut), Tuple{Ptr{GLTriangleFace},Cint},
                  (Ptr{Ptr{Float64}}, Ptr{UInt32}, UInt32), polygon, lengths, len)
    return unsafe_wrap(Vector{GLTriangleFace}, array[1], array[2]; own=true)
end

function earcut_triangulate(polygon::Vector{Vector{Point{2,Float32}}})
    lengths = map(x -> UInt32(length(x)), polygon)
    len = UInt32(length(lengths))
    array = ccall((:u32_triangulate_f32, libearcut), Tuple{Ptr{GLTriangleFace},Cint},
                  (Ptr{Ptr{Float32}}, Ptr{UInt32}, UInt32), polygon, lengths, len)
    return unsafe_wrap(Vector{GLTriangleFace}, array[1], array[2]; own=true)
end

function earcut_triangulate(polygon::Vector{Vector{Point{2,Int64}}})
    lengths = map(x -> UInt32(length(x)), polygon)
    len = UInt32(length(lengths))
    array = ccall((:u32_triangulate_i64, libearcut), Tuple{Ptr{GLTriangleFace},Cint},
                  (Ptr{Ptr{Int64}}, Ptr{UInt32}, UInt32), polygon, lengths, len)
    return unsafe_wrap(Vector{GLTriangleFace}, array[1], array[2]; own=true)
end

function earcut_triangulate(polygon::Vector{Vector{Point{2,Int32}}})
    lengths = map(x -> UInt32(length(x)), polygon)
    len = UInt32(length(lengths))
    array = ccall((:u32_triangulate_i32, libearcut), Tuple{Ptr{GLTriangleFace},Cint},
                  (Ptr{Ptr{Int32}}, Ptr{UInt32}, UInt32), polygon, lengths, len)
    return unsafe_wrap(Vector{GLTriangleFace}, array[1], array[2]; own=true)
end

best_earcut_eltype(x) = Float64
best_earcut_eltype(::Type{Float64}) = Float64
best_earcut_eltype(::Type{<:AbstractFloat}) = Float32
best_earcut_eltype(::Type{Int64}) = Int64
best_earcut_eltype(::Type{Int32}) = Int32
best_earcut_eltype(::Type{<:Integer}) = Int64

function faces(polygon::Polygon{Dim,T}) where {Dim,T}
    PT = Point{Dim,best_earcut_eltype(T)}
    points = [decompose(PT, polygon.exterior)]
    foreach(x -> push!(points, decompose(PT, x)), polygon.interiors)
    return earcut_triangulate(points)
end
