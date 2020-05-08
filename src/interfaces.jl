"""
    coordinates(geometry)
Returns the edges/vertices/coordinates of a geometry. Is allowed to return lazy iterators!
Use `decompose(ConcretePointType, geometry)` to get `Vector{ConcretePointType}` with
`ConcretePointType` to be something like `Point{3, Float32}`.
"""
function coordinates(points::AbstractVector{<:AbstractPoint})
    return points
end

"""
    faces(geometry)
Returns the face connections of a geometry. Is allowed to return lazy iterators!
Use `decompose(ConcreteFaceType, geometry)` to get `Vector{ConcreteFaceType}` with
`ConcreteFaceType` to be something like `TriangleFace{Int}`.
"""
function faces(f::AbstractVector{<:AbstractFace})
    return f
end

coordinates(x, ::Nothing) = coordinates(x)
faces(x, ::Nothing) = faces(x)
decompose(x, ::Nothing) = decompose(x)
decompose(t::Type{<:Point}, x, ::Nothing) = decompose(t, x)
decompose(t::Type{<:NgonFace}, x, ::Nothing) = decompose(t, x)
decompose(t::Vec, x, ::Nothing) = decompose(t, x)
decompose(t::UV, x, ::Nothing) = decompose(t, x)
decompose(t::UVW, x, ::Nothing) = decompose(t, x)
decompose(t::Normal, x, ::Nothing) = decompose(t, x)
