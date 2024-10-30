"""
    coordinates(geometry)

Returns the positions/coordinates of a geometry. 

This is allowed to return lazy iterators. Use `decompose(ConcretePointType, geometry)` 
to get a `Vector{ConcretePointType}` with `ConcretePointType` being something like
`Point3f`.
"""
function coordinates(points::AbstractVector{<:Point})
    return points
end

"""
    faces(geometry)

Returns the faces of a geometry. 

This is allowed to return lazy iterators. Use `decompose(ConcreteFaceType, geometry)` 
to get a `Vector{ConcreteFaceType}` with `ConcreteFaceType` being something like `GLTriangleFace`.
"""
function faces(f::AbstractVector{<:AbstractFace})
    return f
end

"""
    normals(primitive)

Returns the normals of a geometry. 

This is allowed to return lazy iterators. Use `decompose_normals(ConcreteVecType, geometry)` 
to get a `Vector{ConcreteVecType}` with `ConcreteVecType` being something like `Vec3f`.
"""
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

"""
    texturecoordinates(primitive)

Returns the texturecoordinates of a geometry. 

This is allowed to return lazy iterators. Use `decompose_uv(ConcreteVecType, geometry)` 
(or `decompose_uvw`) to get a `Vector{ConcreteVecType}` with `ConcreteVecType` being 
something like `Vec2f`.
"""
texturecoordinates(primitive, nvertices=nothing) = nothing

"""
    Tessellation(primitive, nvertices)

When generating a mesh from an abstract geometry, we can typically generate it
at different levels of detail, i.e. with different amounts of vertices. The 
`Tessellation` wrapper allows you to specify this level of detail. When generating
a mesh from a tessellated geometry, the added information will be passed to 
`coordinates`, `faces`, etc.

```julia
sphere = Sphere(Point3f(0), 1)
m1 = mesh(sphere) # uses a default value for tessellation
m2 = mesh(Tessellation(sphere, 64)) # uses 64 for tessellation
length(coordinates(m1)) != length(coordinates(m2))
```

For grid based tessellation, you can also use a tuple:

```julia
rect = Rect2(0, 0, 1, 1)
Tessellation(rect, (5, 5))
```
"""
struct Tessellation{Dim,T,Primitive,NGrid} <: AbstractGeometry{Dim, T}
    primitive::Primitive
    nvertices::NTuple{NGrid,Int}
end

Base.@deprecate_binding Tesselation Tessellation

function Tessellation(primitive::GeometryPrimitive{Dim,T},
                     nvertices::NTuple{N,<:Integer}) where {Dim,T,N}
    return Tessellation{Dim,T,typeof(primitive),N}(primitive, Int.(nvertices))
end

Tessellation(primitive, nvertices::Integer) = Tessellation(primitive, (nvertices,))

# This is a bit lazy, I guess we should just refactor these methods
# to directly work on Tessellation - but this way it's backward compatible and less
# refactor work :D
nvertices(tessellation::Tessellation) = tessellation.nvertices
nvertices(tessellation::Tessellation{T,N,P,1}) where {T,N,P} = tessellation.nvertices[1]

function coordinates(tessellation::Tessellation)
    return coordinates(tessellation.primitive, nvertices(tessellation))
end
faces(tessellation::Tessellation) = faces(tessellation.primitive, nvertices(tessellation))
normals(tessellation::Tessellation) = normals(tessellation.primitive, nvertices(tessellation))
function texturecoordinates(tessellation::Tessellation)
    return texturecoordinates(tessellation.primitive, nvertices(tessellation))
end

## Decompose methods
# Dispatch type to make `decompose(UV{Vec2f}, primitive)` work
# and to pass through tessellation information

struct UV{T} end
UV(::Type{T}) where {T} = UV{T}()
UV() = UV(Vec2f)
struct UVW{T} end
UVW(::Type{T}) where {T} = UVW{T}()
UVW() = UVW(Vec3f)
struct Normal{T} end
Normal(::Type{T}) where {T} = Normal{T}()
Normal() = Normal(Vec3f)

"""
    decompose(::Type{TargetType}, primitive)
    decompose(::Type{TargetType}, data::AbstractVector)

Dependent on the given type, extracts data from the primtive and converts its 
eltype to `TargetType`.

Possible `TargetType`s:
- `<: Point` extracts and converts positions (calling `coordinates()`)
- `<: AbstractFace` extracts and converts faces (calling `faces()`)
- `<: Normal{<: Vec}` extracts and converts normals, potentially generating them (calling `normals()`)
- `<: UV{<: Vec}` extracts and converts 2D texture coordinates, potentially generating them (calling `texturecoordinates()`)
- `<: UVW{<: Vec}` extracts and converts 3D texture coordinates, potentially generating them (calling `texturecoordinates()`)
"""
function decompose(::Type{F}, primitive::AbstractGeometry) where {F<:AbstractFace}
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
  
function decompose(::Type{F}, f::AbstractVector) where {F<:AbstractFace}
    fs = faces(f)
    isnothing(fs) && error("No faces defined for $(typeof(f))")
    return collect_with_eltype(F, fs)
end

# TODO: Should this be a completely different function?
function decompose(::Type{F}, f::AbstractVector, views::Vector{UnitRange{Int}}) where {F<:AbstractFace}
    fs = faces(f)
    isnothing(fs) && error("No faces defined for $(typeof(f))")
    if isempty(views)
        return collect_with_eltype(F, fs), views
    else
        output = F[]
        new_views = sizehint!(UnitRange{Int}[], length(views))
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

function decompose(::Type{Point{Dim}}, primitive::AbstractGeometry{Dim,T}) where {Dim,T}
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
