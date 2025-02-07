"""
    TupleView{T, N, Skip, A}

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
struct TupleView{T,N,Skip,A} <: AbstractVector{T}
    data::A
    connect::Bool
end

coordinates(tw::TupleView) = coordinates(tw.data)

function Base.size(x::TupleView{T,N,M}) where {T,N,M}
    nitems = length(x.data) ÷ (N - (N - M))
    nitems = nitems - max(N - M, 0)
    return (nitems + x.connect,) # plus one item if we connect
end

function Base.getindex(x::TupleView{T,N,M}, index::Integer) where {T,N,M}
    return ntuple(i -> x.data[mod1(((index - 1) * M) + i, length(x.data))], N)
end

function TupleView{N}(x::AbstractVector; connect=false) where {N}
    return TupleView{N,N}(x, connect=connect)
end

function TupleView{N,M}(x::AbstractVector{T}; connect=false) where {T,N,M}
    return TupleView{NTuple{N,T},N,M,typeof(x)}(x, connect)
end

"""
    connect(points::AbstractVector{<: Point}, P::Type{<: Polytope{N}}, skip::Int = N)

Creates a view that connects a number of points to a Polytope `P`.
Between each polytope, `skip` elements are skipped until the next starts.
Example:
```julia
x = connect(Point[(1, 2), (3, 4), (5, 6), (7, 8)], Line, 2)
x == [Line(Point(1, 2), Point(3, 4)), Line(Point(5, 6), Point(7, 8))]
"""
function connect(points::AbstractVector{Point},
                         P::Type{<:Polytope{N,T} where {N,T}},
                         skip::Int=length(P)) where {Point <: Point}
    return map(Polytope(P, Point), TupleView{length(P),skip}(points))
end

function connect(points::AbstractVector{T}, ::Type{<:Point{N}}, skip::Int=N) where {T <: Real,N}
    return connect(points, Point{N, T}, skip)
end
function connect(points::AbstractVector{_T}, ::Type{<:Point{N, T}}, skip::Int=N) where {T <: Real, N, _T <: Real}
    return map(Point{N,T}, TupleView{N,skip}(points))
end

function connect(indices::AbstractVector{T}, P::Type{<:AbstractFace{N}},
                         skip::Int=N) where {T <: Integer, N}
    return connect(indices, Face(P, T), skip)
end
function connect(indices::AbstractVector{_T}, P::Type{<:AbstractFace{N, T}},
                         skip::Int=N) where {T <: Integer, N, _T <: Integer}
    return collect(reinterpret(Face(P, T), TupleView{N, skip}(indices)))
end

function connect(points::AbstractMatrix{T}, ::Type{<:Point{N}}) where {T <: Real, N}
    return connect(points, Point{N, T})
end
function connect(points::AbstractMatrix{_T}, P::Type{Point{N, T}}) where {T <: Real, N, _T <: Real}
    return if size(points, 1) === N
        return reinterpret(Point{N,T}, points)
    elseif size(points, 2) === N
        seglen = size(points, 1)
        columns = ntuple(N) do i
            return view(points, ((i - 1) * seglen + 1):(i * seglen))
        end
        return P.(columns...)
    else
        error("Dim 1 or 2 must be equal to the point dimension!")
    end
end

function connect(elements::AbstractVector, faces::AbstractVector{<: AbstractFace})
    return [elements[i] for f in faces for i in f]
end
