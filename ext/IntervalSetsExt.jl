module IntervalSetsExt

using IntervalSets
using GeometryBasics

GeometryBasics.HyperRectangle(ints::Vararg{ClosedInterval, N}) where {N} = HyperRectangle{N}(ints...)
GeometryBasics.HyperRectangle{N}(ints::Vararg{ClosedInterval, N}) where {N} = HyperRectangle{N}(
    Vec(leftendpoint.(ints)),
    Vec(rightendpoint.(ints) .- leftendpoint.(ints))
)
GeometryBasics.HyperRectangle{N,T}(ints::Vararg{ClosedInterval, N}) where {N,T} = HyperRectangle{N,T}(
    Vec(leftendpoint.(ints)),
    Vec(rightendpoint.(ints) .- leftendpoint.(ints))
)

end
