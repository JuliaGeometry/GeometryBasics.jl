"""
Abstract Geometry in R{Dim} with Number type T
"""
abstract type AbstractGeometry{Dim,T<:Number} end
abstract type GeometryPrimitive{Dim,T} <: AbstractGeometry{Dim,T} end
Base.ndims(x::AbstractGeometry{Dim}) where {Dim} = Dim

"""
Geometry made of N connected points. Connected as one flat geometry, it makes a Ngon / Polygon.
Connected as volume it will be a Simplex / Tri / Cube.
Note That `Polytope{N} where N == 3` denotes a Triangle both as a Simplex or Ngon.
"""
abstract type Polytope{Dim,T} <: AbstractGeometry{Dim,T} end
abstract type AbstractPolygon{Dim,T} <: Polytope{Dim,T} end

abstract type AbstractFace{N,T} <: StaticVector{N,T} end
abstract type AbstractSimplexFace{N,T} <: AbstractFace{N,T} end
abstract type AbstractNgonFace{N,T} <: AbstractFace{N,T} end

abstract type AbstractSimplex{Dim,N,T} <: StaticVector{Dim,T} end

@fixed_vector SimplexFace AbstractSimplexFace

const TetrahedronFace{T} = SimplexFace{4,T}
Face(::Type{<:SimplexFace{N}}, ::Type{T}) where {N,T} = SimplexFace{N,T}

@fixed_vector NgonFace AbstractNgonFace

const LineFace{T} = NgonFace{2,T}
const TriangleFace{T} = NgonFace{3,T}
const QuadFace{T} = NgonFace{4,T}

function Base.show(io::IO, x::TriangleFace{T}) where {T}
    return print(io, "TriangleFace(", join(x, ", "), ")")
end

Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F

@propagate_inbounds Base.getindex(x::Polytope, i::Integer) = coordinates(x)[i]
@propagate_inbounds Base.iterate(x::Polytope) = iterate(coordinates(x))
@propagate_inbounds Base.iterate(x::Polytope, i) = iterate(coordinates(x), i)

"""
Fixed Size Polygon, e.g.
- N 1-2 : Illegal!
- N = 3 : Triangle
- N = 4 : Quadrilateral (or Quad, Or tetragon)
- N = 5 : Pentagon
- ...
"""
struct Ngon{Dim, T<:Real, N} <: AbstractPolygon{Dim,T}
    points::NTuple{N, Point{Dim, T}}
end

const NNgon{N} = Ngon{Dim,T,N} where {Dim,T}

function (::Type{<: NNgon{N}})(points::Vararg{Point{Dim,T}, N}) where {N,Dim,T}
    return Ngon{Dim,T,N}(points)
end
Base.show(io::IO, x::NNgon{N}) where {N} = print(io, "Ngon{$N}(", join(x, ", "), ")")

# Interfaces
coordinates(x::Ngon) = x.points
# Base Array interface
Base.length(::Type{<:NNgon{N}}) where {N} = N
Base.length(::NNgon{N}) where {N} = N

"""
The Ngon Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{Point{Dim,T}},
                  ::Type{<:AbstractNgonFace{N,IT}}) where {N,Dim,T,IT}
    return Ngon{Dim,T,N,P}
end

"""
The fully concrete Ngon type, when constructed from a point type!
"""
function Polytope(::Type{<:NNgon{N}}, P::Type{Point{NDim,T}}) where {N,NDim,T}
    return Ngon{NDim,T,N}
end

const Line{Dim,T} = Ngon{Dim,T,2}

# Simplex{D, T, 3} & Ngon{D, T, 3} are both representing a triangle.
# Since Ngon is supposed to be flat and a triangle is flat, lets prefer Ngon
# for triangle:
const Triangle{Dim,T} = Ngon{Dim,T,3}
const Triangle3d{T} = Triangle{3,T}

Base.show(io::IO, x::Triangle) = print(io, "Triangle(", join(x, ", "), ")")

const Quadrilateral{Dim,T} = Ngon{Dim,T,4}

Base.show(io::IO, x::Quadrilateral) = print(io, "Quad(", join(x, ", "), ")")

function coordinates(lines::AbstractArray{Line{Dim,T}}) where {Dim,T}
    result = Point{Dim, T}[]
    for line in lines
        append!(result, coordinates(line))
    end
    return result
end

"""
A `Simplex` is a generalization of an N-dimensional tetrahedra and can be thought
of as a minimal convex set containing the specified points.

* A 0-simplex is a point.
* A 1-simplex is a line segment.
* A 2-simplex is a triangle.
* A 3-simplex is a tetrahedron.

Note that this datatype is offset by one compared to the traditional
mathematical terminology. So a one-simplex is represented as `Simplex{2,T}`.
This is for a simpler implementation.

It applies to infinite dimensions. The structure of this type is designed
to allow embedding in higher-order spaces by parameterizing on `T`.
"""
struct Simplex{Dim,T<:Real,N} <: Polytope{Dim,T}
    points::NTuple{N,Point{Dim,T}}
end

const NSimplex{N} = Simplex{Dim,T,N} where {Dim,T}
const Tetrahedron{T} = Simplex{3,T,4}

Base.show(io::IO, x::Tetrahedron) = print(io, "Tetrahedron(", join(x, ", "), ")")

coordinates(x::Simplex) = x.points

function (::Type{<:NSimplex{N}})(points::Vararg{Point{Dim,T},N}) where {Dim,T,N}
    return Simplex{Dim,T,N}(points)
end

# Base Array interface
Base.length(::Type{<:NSimplex{N}}) where {N} = N
Base.length(::NSimplex{N}) where {N} = N

"""
The Simplex Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{Point{Dim,T}}, ::Type{<:AbstractSimplexFace{N}}) where {N,Dim,T}
    return Simplex{Dim,T,N,P}
end

"""
The fully concrete Simplex type, when constructed from a point type!
"""
function Polytope(::Type{<:NSimplex{N}}, P::Type{Point{NDim,T}}) where {N,NDim,T}
    return Simplex{NDim,T,N,P}
end
Base.show(io::IO, x::Line) = print(io, "Line(", x[1], " => ", x[2], ")")

"""
    LineString(points::AbstractVector{<:Point})

A LineString is a geometry of connected line segments
"""
struct LineString{Dim,T<:Real,V<:AbstractVector{Line{Dim,T}}}
    points::V
end

coordinates(x::LineString) = coordinates(x.points)

Base.copy(x::LineString) = LineString(copy(x.points))
Base.size(x::LineString) = size(getfield(x, :points))
Base.getindex(x::LineString, i) = getindex(getfield(x, :points), i)

function LineString(points::AbstractVector{Line{Dim,T}}) where {Dim,T}
    return LineString{Dim,T,typeof(points)}(points)
end

"""
    LineString(points::AbstractVector{<: Point}, skip = 1)

Creates a LineString from a vector of points.
With `skip == 1`, the default, it will connect the line like this:
```julia
points = Point[a, b, c, d]
linestring = LineString(points)
@assert linestring == LineString([a => b, b => c, c => d])
```
"""
function LineString(points::AbstractVector{<:Point}, skip=1)
    return LineString(connect(points, Line, skip))
end

function LineString(points::AbstractVector{<:Pair{Point{N,T},Point{N,T}}}) where {N, T}
    return LineString(reinterpret(Line{N,T}, points))
end

function LineString(points::AbstractVector{<:Point},
                    faces::AbstractVector{<:LineFace})
    return LineString(connect(points, faces))
end

"""
    LineString(points::AbstractVector{<: Point}, indices::AbstractVector{<: Integer}, skip = 1)

Creates a LineString from a vector of points and an index list.
With `skip == 1`, the default, it will connect the line like this:


```julia
points = Point[a, b, c, d]; faces = [1, 2, 3, 4]
linestring = LineString(points, faces)
@assert linestring == LineString([a => b, b => c, c => d])
```
To make a segmented line, set skip to 2
```julia
points = Point[a, b, c, d]; faces = [1, 2, 3, 4]
linestring = LineString(points, faces, 2)
@assert linestring == LineString([a => b, c => d])
```
"""
function LineString(points::AbstractVector{<:Point},
                    indices::AbstractVector{<:Integer}, skip=1)
    faces = connect(indices, LineFace, skip)
    return LineString(points, faces)
end

"""
    Polygon(exterior::AbstractVector{<:Point})
    Polygon(exterior::AbstractVector{<:Point}, interiors::Vector{<:AbstractVector{<:Point}})

"""
struct Polygon{Dim,T<:Real,L<:AbstractVector{Point{Dim,T}},
               V<:AbstractVector{L}} <: AbstractPolygon{Dim,T}
    exterior::L
    interiors::V
end

Base.copy(x::Polygon) = Polygon(copy(x.exterior), copy(x.interiors))

function Base.:(==)(a::Polygon, b::Polygon)
    return (a.exterior == b.exterior) && (a.interiors == b.interiors)
end

function Polygon(
        exterior::AbstractVector{Point{Dim,T}},
        interiors::AbstractVector{AbstractVector{Point{Dim,T}}}) where {Dim, T}
    return Polygon{Dim,T,typeof(exterior),typeof(interiors)}(exterior, interiors)
end

Polygon(exterior::AbstractVector{Point{N, T}}) where {N, T} = Polygon(exterior, Vector{Point{N, T}}[])

function Polygon(exterior::AbstractVector{Point{Dim, T}}, faces::AbstractVector{<:Integer},
                 skip::Int=1) where {Dim,T}
    return Polygon(LineString(exterior, faces, skip))
end

function Polygon(exterior::AbstractVector{Point{Dim,T}},
                 faces::AbstractVector{<:LineFace}) where {Dim, T}
    return Polygon(LineString(exterior, faces))
end

function Polygon(exterior::AbstractVector{Point{Dim, T}},
                 interior::AbstractVector{<:AbstractVector{Point{Dim, T}}}) where {Dim,T}
    ext = LineString(exterior)
    # We need to take extra care for empty interiors, since
    # if we just map over it, it won't infer the element type correctly!
    int = typeof(ext)[]
    foreach(x -> push!(int, LineString(x)), interior)
    return Polygon(ext, int)
end

function coordinates(polygon::Polygon{N,T}) where {N,T}
    exterior = coordinates(polygon.exterior)
    if isempty(polygon.interiors)
        return exterior
    else
        result = Point{N, T}[]
        append!(result, exterior)
        foreach(x -> append!(result, coordinates(x)), polygon.interiors)
        return result
    end
end

"""
    MultiPolygon(polygons::AbstractPolygon)
"""
struct MultiPolygon{Dim,T<:Real,A<:AbstractVector{<:AbstractPolygon{Dim,T}}}
    polygons::A
end

function MultiPolygon(polygons::AbstractVector{<: AbstractPolygon{Dim,T}}) where {Dim,T}
    return MultiPolygon(polygons)
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)

struct MultiLineString{Dim,T<:Real,A<:AbstractVector{LineString{Dim,T}}}
    linestrings::A
end

function MultiLineString(linestrings::AbstractVector{L}) where {L<:AbstractVector{Line{Dim,T}}} where {Dim,T}
    return MultiLineString(linestrings)
end

Base.getindex(ms::MultiLineString, i) = ms.linestrings[i]
Base.size(ms::MultiLineString) = size(ms.linestrings)

"""
    MultiPoint(points::AbstractVector{<: Point})

A collection of points
"""
struct MultiPoint{Dim,T<:Real,A<:AbstractVector{Point{Dim,T}}}
    points::A
end

function MultiPoint(points::AbstractVector{Point{Dim,T}}) where {Dim,T}
    return MultiPoint(points)
end

Base.getindex(mpt::MultiPoint, i) = mpt.points[i]
Base.size(mpt::MultiPoint) = size(mpt.points)

"""
    AbstractMesh

An abstract mesh is a collection of Polytope elements (Simplices / Ngons).
The connections are defined via faces(mesh), the coordinates of the elements are returned by
coordinates(mesh). Arbitrary meta information can be attached per point or per face
"""
abstract type AbstractMesh{Element} end

"""
    Mesh <: AbstractVector{Element}
The concrete AbstractMesh implementation.
"""
struct Mesh{Dim,T<:Number, Element, V<:AbstractVector{Point{Dim, T}}, C} <:
       AbstractMesh{Element}
    vertices::V
    connectivity::C
end

function Base.summary(io::IO, ::Mesh{Dim,T,Element}) where {Dim,T,Element}
    print(io, "Mesh{$Dim, $T, ")
    summary(io, Element)
    return print(io, "}")
end

function Mesh(elements::AbstractVector{<:Polytope{Dim,T}}) where {Dim,T}
    return Mesh{Dim,T,eltype(elements),typeof(elements)}(elements, nothing)
end

function Mesh(points::AbstractVector{<:Point},
              faces::AbstractVector{<:AbstractFace})
    return Mesh(points, faces)
end

function Mesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(connect(points, connect(faces, facetype, skip)))
end
