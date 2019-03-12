"""
TupleView, groups elements of an array as tuples.
N is the dimension of the tuple, M tells how many elements to skip to the next tuple.
By default TupleView{N} defaults to skip N items.
# a few examples:

```julia

x = [1, 2, 3, 4, 5, 6]
TupleView{2, 1}(x):
> [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

TupleView{2}(x):
> [(1, 2), (3, 4), (5, 6)]

TupleView{2, 3}(x):
> [(1, 2), (4, 5)]

TupleView{3, 1}(x):
> [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6)]
```

TupleView can be used together with reinterpret:

```julia
x = [1, 2, 3, 4]
y = reinterpret(Point{2, Int}, TupleView{2, 1}(x))
> [Point(1, 2), Point(2, 3), Point(3, 4)]
```

"""
struct TupleView{T, N, Skip, A} <: AbstractVector{T}
    data::A
    connect::Bool
end
Base.show(io::IO, ::Type{<: TupleView{T, N, Skip}}) where {T, N, Skip} = print(io, "TupleView{$T, $Skip}")

function Base.size(x::TupleView{T, N, M}) where {T, N, M}
    nitems = length(x.data) ÷ (N - (N - M))
    nitems = nitems - max(N - M, 0)
    (nitems + x.connect,) # plus one item if we connect
end

function Base.getindex(x::TupleView{T, N, M}, index::Integer) where {T, N, M}
    return ntuple(i-> x.data[mod1(((index - 1) * M) + i, length(x.data))], N)
end

TupleView{N}(x::AbstractVector; connect = false) where N = TupleView{N, N}(x, connect = connect)

function TupleView{N, M}(x::AbstractVector{T}; connect = false) where {T, N, M}
    TupleView{NTuple{N, T}, N, M, typeof(x)}(x, connect)
end


@inline connected_line(points::AbstractVector{<: AbstractPoint{N}}, skip = N) where N = connect(points, Line, skip)

"""
    connect(points::AbstractVector{<: AbstractPoint}, P::Type{<: Polytype{N}}, skip::Int = N)

Creates a view that connects a number of points to a Polytope `P`.
Between each polytope, `skip` elements are skipped untill the next starts.
Example:
```julia
x = connect(Point[(1, 2), (3, 4), (5, 6), (7, 8)], Line, 2)
x == [Line(Point(1, 2), Point(3, 4)), Line(Point(5, 6), Point(7, 8))]
"""
@inline function connect(points::AbstractVector{Point}, P::Type{<: Polytope{N, T} where {N, T}}, skip::Int = length(P)) where Point <: AbstractPoint
    reinterpret(Polytope(P, Point), TupleView{length(P), skip}(points))
end

@inline function connect(points::AbstractVector{T}, P::Type{<: Point{N}}, skip::Int = N) where {T <: Real, N}
    reinterpret(Point{N, T}, TupleView{N, skip}(points))
end

@inline function connect(points::AbstractVector{T}, P::Type{<: AbstractFace{N}}, skip::Int = N) where {T <: Real, N}
    reinterpret(Face(P, T), TupleView{N, skip}(points))
end

"""
FaceView enables to link one array of points via a face array, to generate one
abstract array of elements.
E.g., this becomes possible:
```
x = FaceView(rand(Point3f0, 10), TriangleFace[(1, 2, 3), (2, 4, 5), ...])
x[1] isa Triangle == true
x isa AbstractVector{<: Triangle} == true
# This means we can use it as a mesh:
Mesh(x) # should just work!
Can also be used for Points:

linestring = FaceView(points, LineFace[...])
Polygon(linestring)
```
"""
struct FaceView{
        Element,
        Point <: AbstractPoint,
        Face <: AbstractFace,
        P <: AbstractVector{Point},
        F <: AbstractVector{Face}
    } <: AbstractVector{Element}

    points::P
    faces::F
end

Base.size(x::FaceView) = size(x.faces)

Base.show(io::IO, x::Type{<: FaceView{Element}}) where Element = print(io, "FaceView{", Element, "}")


@propagate_inbounds function Base.getindex(x::FaceView{Element}, i) where Element
    Element(map(idx-> x.points[idx], x.faces[i]))
end

@propagate_inbounds function Base.setindex!(x::FaceView{Element}, element::Element, i) where Element
    face = x.faces[i]
    for (i, f) in enumerate(face) # TODO unroll!?
        x.points[face[i]] = element[i]
    end
    return element
end

function connect(points::AbstractVector{P}, faces::AbstractVector{F}) where {P <: AbstractPoint, F <: AbstractFace}
    FaceView{Polytope(P, F), P, F, typeof(points), typeof(faces)}(points, faces)
end



const FaceMesh{Dim, T, Element} = Mesh{Dim, T, Element, <: FaceView{Element}}

function coordinates(mesh::FaceMesh)
    mesh.simplices.points
end
function faces(mesh::FaceMesh)
    mesh.simplices.faces
end
