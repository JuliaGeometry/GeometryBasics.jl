"""
    coordinates(geometry)
Returns the edges/vertices/coordinates of a geometry. Is allowed to return lazy iterators!
Use `decompose(ConcretePointType, geometry)` to get `Vector{ConcretePointType}` with
`ConcretePointType` to be something like `Point{3, Float32}`.
"""
function coordinates(points::AbstractVector{<:Point})
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
sphere = Sphere(Point3f(0), 1)
m1 = mesh(sphere) # uses a default value for tesselation
m2 = mesh(Tesselation(sphere, 64)) # uses 64 for tesselation
length(coordinates(m1)) != length(coordinates(m2))
```
For grid based tesselation, you can also use a tuple:
```julia
rect = Rect2(0, 0, 1, 1)
Tesselation(rect, (5, 5))
"""
struct Tesselation{Dim,T,Primitive,NGrid} <: AbstractGeometry{Dim, T}
    primitive::Primitive
    nvertices::NTuple{NGrid,Int}
end

function Tesselation(primitive::GeometryPrimitive{Dim,T},
                     nvertices::NTuple{N,<:Integer}) where {Dim,T,N}
    return Tesselation{Dim,T,typeof(primitive),N}(primitive, Int.(nvertices))
end

Tesselation(primitive, nvertices::Integer) = Tesselation(primitive, (nvertices,))

# This is a bit lazy, I guess we should just refactor these methods
# to directly work on Tesselation - but this way it's backward compatible and less
# refactor work :D
nvertices(tesselation::Tesselation) = tesselation.nvertices
nvertices(tesselation::Tesselation{T,N,P,1}) where {T,N,P} = tesselation.nvertices[1]

function coordinates(tesselation::Tesselation)
    return coordinates(tesselation.primitive, nvertices(tesselation))
end
faces(tesselation::Tesselation) = faces(tesselation.primitive, nvertices(tesselation))
normals(tesselation::Tesselation) = normals(tesselation.primitive, nvertices(tesselation))
function texturecoordinates(tesselation::Tesselation)
    return texturecoordinates(tesselation.primitive, nvertices(tesselation))
end

## Decompose methods
# Dispatch type to make `decompose(UV{Vec2f}, primitive)` work
# and to pass through tesselation information

struct UV{T} end
UV(::Type{T}) where {T} = UV{T}()
UV() = UV(Vec2f)
struct UVW{T} end
UVW(::Type{T}) where {T} = UVW{T}()
UVW() = UVW(Vec3f)
struct Normal{T} end
Normal(::Type{T}) where {T} = Normal{T}()
Normal() = Normal(Vec3f)

function decompose(::Type{F}, primitive::AbstractGeometry) where {F<:AbstractVertexFace}
    f = faces(primitive)
    if isnothing(f)
        if ndims(primitive) == 2
            # if 2d, we can fallback to Polygon triangulation
            return decompose(F, Polygon(decompose(Point, primitive)))
        else
            return nothing
        end
    end
    return decompose(F, f)
end
  
function decompose(::Type{F}, f::AbstractVector) where {F<:AbstractVertexFace}
    fs = faces(f)
    isnothing(fs) && error("No faces defined for $(typeof(f))")
    return collect_with_eltype(F, fs)
end

# TODO: Should this be a completely different function?
function decompose(::Type{F}, f::AbstractVector, views::Vector{UnitRange{Int}}) where {F<:AbstractVertexFace}
    fs = faces(f)
    isnothing(fs) && error("No faces defined for $(typeof(f))")
    if isempty(views)
        return collect_with_eltype(F, fs), views
    else
        output = F[]
        new_views = UnitRange{Int}[]
        for range in views
            start = length(output) + 1
            collect_with_eltype!(output, view(fs, range))
            push!(new_views, start:length(output))
        end
        return output, new_views
    end
end

function decompose(::Type{P}, primitive) where {P<:Point}
    return collect_with_eltype(P, coordinates(primitive))
end

function decompose(::Type{Point}, primitive::AbstractGeometry{Dim,T}) where {Dim,T}
    return collect_with_eltype(Point{Dim,T}, coordinates(primitive))
end

function decompose(::Type{PointT{T}}, primitive::AbstractGeometry{Dim}) where {Dim, T}
    return collect_with_eltype(Point{Dim,T}, coordinates(primitive))
end

function decompose(::Type{T}, primitive) where {T}
    return collect_with_eltype(T, primitive)
end

decompose_uv(primitive) = decompose(UV(), primitive)
decompose_uvw(primitive) = decompose(UVW(), primitive)
decompose_normals(primitive) = decompose(Normal(), primitive)

decompose_uv(T, primitive) = decompose(UV(T), primitive)
decompose_uvw(T, primitive) = decompose(UVW(T), primitive)
decompose_normals(T, primitive) = decompose(Normal(T), primitive)

function decompose(::Normal{T}, primitive) where {T}
    _n = normals(primitive)
    if isnothing(_n)
        # For 3D primitives, we can calculate the normals from the vertices + faces
        if ndims(primitive) == 3
            n = normals(decompose(Point, primitive), faces(primitive), T)
        else
            points = decompose(Point, primitive)
            n = [T(0, 0, 1) for p in points]
        end
    else
        n = _n
    end
    return collect_with_eltype(T, n)
end

function decompose(UVT::Union{UV{T},UVW{T}}, primitive) where {T}
    # This is the fallback for texture coordinates if a primitive doesn't overload them
    # We just take the positions and normalize them
    uv = texturecoordinates(primitive)
    if uv === nothing
        # If the primitive doesn't even have coordinates, we're out of options and return
        # nothing, indicating that texturecoordinates aren't implemented
        positions = decompose(Point, primitive)
        isnothing(positions) && return nothing
        # We should generate 3D uvw's if uv's are requested
        # TODO: we should probably enforce that UV has a 2D type and UVW a 3D type
        (length(T) != length(eltype(positions))) && return nothing

        # Let this overlord do the work
        return decompose(UVT, positions)
    end
    return collect_with_eltype(T, uv)
end

function decompose(::Union{UV{T},UVW{T}}, positions::AbstractVector{<:VecTypes}) where {T}
    N = length(T)
    positions_nd = decompose(Point{N,eltype(T)}, positions)
    bb = Rect(positions_nd) # Make sure we get this as points
    mini, w = minimum(bb), widths(bb)
    return map(positions_nd) do p
        return T((p .- mini) ./ w)
    end
end
