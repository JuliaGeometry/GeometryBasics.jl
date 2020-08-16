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

function normals(primitive, nvertices=nothing; kw...)
    # doesn't have any specific algorithm to generate normals
    # so will be generated from faces + positions
    # which we indicate by returning nothing!
    # Overload normals(primitive::YourPrimitive), to calcalute the normals
    # differently
    return nothing
end

function faces(primitive, nvertices=nothing; kw...)
    # doesn't have any specific algorithm to generate faces
    # so will try to triangulate the coordinates!
    return nothing
end

texturecoordinates(primitive, nvertices=nothing) = nothing

"""
    Tesselation(primitive, nvertices)
For abstract geometries, when we generate
a mesh from them, we need to decide how fine grained we want to mesh them.
To transport this information to the various decompose methods, you can wrap it
in the Tesselation object e.g. like this:

```julia
sphere = Sphere(Point3f0(0), 1)
m1 = mesh(sphere) # uses a default value for tesselation
m2 = mesh(Tesselation(sphere, 64)) # uses 64 for tesselation
length(coordinates(m1)) != length(coordinates(m2))
```
For grid based tesselation, you can also use a tuple:
```julia
rect = Rect2D(0, 0, 1, 1)
Tesselation(rect, (5, 5))
"""
struct Tesselation{Dim, T, Primitive, NGrid}
    primitive::Primitive
    nvertices::NTuple{NGrid, Int}
end

function Tesselation(primitive::GeometryPrimitive{Dim, T}, nvertices::NTuple{N, <:Integer}) where {Dim, T, N}
    return Tesselation{Dim, T, typeof(primitive), N}(primitive, Int.(nvertices))
end

Tesselation(primitive, nvertices::Integer) = Tesselation(primitive, (nvertices,))

# This is a bit lazy, I guess we should just refactor these methods
# to directly work on Tesselation - but this way it's backward compatible and less
# refactor work :D
nvertices(tesselation::Tesselation) = tesselation.nvertices
nvertices(tesselation::Tesselation{T, N, P, 1}) where {T, N, P} = tesselation.nvertices[1]

coordinates(tesselation::Tesselation) = coordinates(tesselation.primitive, nvertices(tesselation))
faces(tesselation::Tesselation) = faces(tesselation.primitive, nvertices(tesselation))
normals(tesselation::Tesselation) = normals(tesselation.primitive, nvertices(tesselation))
texturecoordinates(tesselation::Tesselation) = texturecoordinates(tesselation.primitive, nvertices(tesselation))

## Decompose methods
# Dispatch type to make `decompose(UV{Vec2f0}, primitive)` work
# and to pass through tesselation information

# Types that can be converted to a mesh via the functions below
const Meshable{Dim, T} = Union{Tesselation{Dim, T}, Mesh{Dim, T}, AbstractPolygon{Dim, T},
                               GeometryPrimitive{Dim, T}, AbstractVector{<: AbstractPoint{Dim, T}}}

struct UV{T} end
UV(::Type{T}) where T = UV{T}()
UV() = UV(Vec2f0)
struct UVW{T} end
UVW(::Type{T}) where T = UVW{T}()
UVW() = UVW(Vec3f0)
struct Normal{T} end
Normal(::Type{T}) where T = Normal{T}()
Normal() = Normal(Vec3f0)

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    f = faces(primitive)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

function decompose(::Type{P}, primitive) where {P<:AbstractPoint}
    return collect_with_eltype(P, metafree(coordinates(primitive)))
end

function decompose(::Type{Point}, primitive::Meshable{Dim, T}) where {Dim, T}
    return collect_with_eltype(Point{Dim, T}, metafree(coordinates(primitive)))
end

function decompose(::Type{Point}, primitive::LineString{Dim, T}) where {Dim, T}
    return collect_with_eltype(Point{Dim, T}, metafree(coordinates(primitive)))
end

function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

decompose_uv(primitive) = decompose(UV(), primitive)
decompose_uvw(primitive) = decompose(UVW(), primitive)
decompose_normals(primitive) = decompose(Normal(), primitive)

function decompose(NT::Normal{T}, primitive) where T
    n = normals(primitive)
    if n === nothing
        return collect_with_eltype(T, normals(coordinates(primitive), faces(primitive)))
    end
    return collect_with_eltype(T, n)
end

function decompose(UVT::Union{UV{T}, UVW{T}}, primitive) where T
    # This is the fallback for texture coordinates if a primitive doesn't overload them
    # We just take the positions and normalize them
    uv = texturecoordinates(primitive)
    if uv === nothing
        # If the primitive doesn't even have coordinates, we're out of options and return
        # nothing, indicating that texturecoordinates aren't implemented
        positions = decompose(Point, primitive)
        positions === nothing && return nothing
        # Let this overlord do the work
        return decompose(UVT, positions)
    end
    return collect_with_eltype(T, uv)
end

function decompose(UVT::Union{UV{T}, UVW{T}}, positions::AbstractVector{<:VecTypes}) where T
    N = length(T)
    positions_nd = decompose(Point{N, eltype(T)}, positions)
    bb = Rect(positions_nd) # Make sure we get this as points
    return map(positions_nd) do p
        return (p .- minimum(bb)) ./ widths(bb)
    end
end

# Stay backward compatible:

decompose(::Type{T}, primitive::Meshable, nvertices) where T = decompose(T, Tesselation(primitive, nvertices))
