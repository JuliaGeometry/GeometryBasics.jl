"""
Abstract Geometry in R{Dim} with Number type T
"""
abstract type AbstractGeometry{Dim, T <: Real} end

"""
Geometry made of N connected points. Connected as one flat geometry, it makes a Ngon / Polygon.
Connected as volume it will be a Simplex / Tri / Cube.
Note That `Polytype{N} where N == 3` denotes a Triangle both as a Simplex or Ngon.
"""
abstract type Polytope{Dim, T} <: AbstractGeometry{Dim, T} end
abstract type AbstractPolygon{Dim, T} <: Polytope{Dim, T} end

abstract type AbstractPoint{Dim, T} <: StaticVector{Dim, T} end
abstract type AbstractFace{Dim, T} <: StaticVector{Dim, T} end
abstract type AbstractSimplexFace{Dim, T} <: AbstractFace{Dim, T} end
abstract type AbstractNgonFace{Dim, T} <: AbstractFace{Dim, T} end

abstract type AbstractSimplex{Dim, N, T} <: StaticVector{Dim, T} end

GeometryTypes.@fixed_vector Point AbstractPoint




"""
Face index, connecting points to form a simplex
"""

GeometryTypes.@fixed_vector SimplexFace AbstractSimplexFace
const LineFace{T} = SimplexFace{2, T}
Base.show(io::IO, x::Type{<: LineFace{T}}) where T = print(io, "LineFace{", T, "}")
Face(::Type{<: SimplexFace{N}}, ::Type{T}) where {N, T} = SimplexFace{N, T}



"""
Face index, connecting points to form an Ngon
"""

GeometryTypes.@fixed_vector NgonFace AbstractNgonFace
const TriangleFace{T} = NgonFace{3, T}
Base.show(io::IO, x::Type{<: TriangleFace{T}}) where T = print(io, "TriangleFace{", T, "}")
Face(::Type{<: NgonFace{N}}, ::Type{T}) where {N, T} = NgonFace{N, T}

@propagate_inbounds Base.getindex(x::Polytope, i::Integer) = coordinates(x)[i]
@propagate_inbounds Base.iterate(x::Polytope) = iterate(coordinates(x))
@propagate_inbounds Base.iterate(x::Polytope, i) = iterate(coordinates(x), i)

"""
Fixed Size Polygon, e.g.
N 1-2 : Illegal!
N = 3 : Triangle
N = 4 : Quadrilateral (or Quad, Or tetragon)
N = 5 : Pentagon
...
"""
struct Ngon{
        Dim, T <: Real,
        N,
        Point <: AbstractPoint{Dim, T}
    } <: AbstractPolygon{Dim, T}

    points::SVector{N, Point}
end
const NNgon{N} = Ngon{Dim, T, N, P} where {Dim, T, P}

function (::Type{<: NNgon{N}})(points::Vararg{P, N}) where {P <: AbstractPoint{Dim, T}, N} where {Dim, T}
    Ngon{Dim, T, N, P}(SVector(points))
end

# Interfaces
coordinates(x::Ngon) = x.points
# Base Array interface
Base.length(::Type{<: NNgon{N}}) where N = N
Base.length(::NNgon{N}) where N = N

"""
The Ngon Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{<: AbstractPoint{Dim, T}}, ::Type{<: AbstractNgonFace{N}}) where {N, Dim, T}
    Ngon{Dim, T, N, P}
end

"""
The fully concrete Ngon type, when constructed from a point type!
"""
function Polytope(::Type{<: NNgon{N}}, P::Type{<: AbstractPoint{NDim, T}}) where {N, NDim, T}
    Ngon{NDim, T, N, P}
end

# Simplex{D, T, 3} & Ngon{D, T, 3} are both representing a triangle.
# Since Ngon is supposed to be flat and a triangle is flat, lets prefer Ngon
# for triangle:
const TriangleP{Dim, T, P <: AbstractPoint{Dim, T}} = Ngon{Dim, T, 3, P}

const Triangle{Dim, T} = TriangleP{Dim, T, Point{Dim, T}}
const Triangle3d{T} = Triangle{3, T}

function Base.show(io::IO, ::Type{TriangleP{D, T, P}}) where {D, T, P}
    print(io, "Tetrahedron{", D, ", ", T, "}")
end
Base.show(io::IO, ::Type{Triangle}) = print(io, "Triangle")
Base.show(io::IO, x::TriangleP) = print(io, "Triangle(", join(x, ", "), ")")

const Quadrilateral{Dim, T} = Ngon{Dim, T, 4, P} where P <: AbstractPoint{Dim, T}

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
struct Simplex{
        Dim, T <: Real,
        N,
        Point <: AbstractPoint{Dim, T},
    } <: Polytope{Dim, T}

    points::SVector{N, Point}
end

const NSimplex{N} = Simplex{Dim, T, N, P} where {Dim, T, P}
const LineP{Dim, T, P <: AbstractPoint{Dim, T}} = Simplex{Dim, T, 2, P}
const Line{Dim, T} = LineP{Dim, T, Point{Dim, T}}
const TetrahedronP{T, P <: AbstractPoint{3, T}} = Simplex{3, T, 4, P}
const Tetrahedron{T} = TetrahedronP{T, Point{3, T}}

function Base.show(io::IO, ::Type{TetrahedronP{T, P}}) where {T, P}
    print(io, "Tetrahedron{", T, "}")
end
Base.show(io::IO, ::Type{Tetrahedron}) = print(io, "Tetrahedron")
Base.show(io::IO, x::TetrahedronP) = print(io, "Tetrahedron(", join(x, ", "), ")")


coordinates(x::Simplex) = x.points
function (::Type{<: NSimplex{N}})(points::Vararg{P, N}) where {P <: AbstractPoint{Dim, T}, N} where {Dim, T}
    Simplex{Dim, T, N, P}(SVector(points))
end

# Base Array interface
Base.length(::Type{<: NSimplex{N}}) where N = N
Base.length(::NSimplex{N}) where N = N

"""
The Simplex Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(P::Type{<: AbstractPoint{Dim, T}}, ::Type{<: AbstractSimplexFace{N}}) where {N, Dim, T}
    Simplex{Dim, T, N, P}
end

"""
The fully concrete Simplex type, when constructed from a point type!
"""
function Polytope(::Type{<: NSimplex{N}}, P::Type{<: AbstractPoint{NDim, T}}) where {N, NDim, T}
    Simplex{NDim, T, N, P}
end

function Base.show(io::IO, ::Type{LineP{D, T, P}}) where {D, T, P}
    print(io, "Line{", D, ", ", T, "}")
end
Base.show(io::IO, ::Type{Line}) = print(io, "Line")
Base.show(io::IO, x::LineP) = print(io, "Line(", x[1], " => ", x[2], ")")

"""
A LineString is a geometry of connected line segments
"""
struct LineString{
        Dim, T <: Real,
        P <: AbstractPoint,
        V <: AbstractVector{<: LineP{Dim, T, P}}
    } <: AbstractVector{LineP{Dim, T}}
    points::V
end
coordinates(x::LineString) = x.points

Base.size(x::LineString) = size(coordinates(x))
Base.getindex(x::LineString, i) = getindex(coordinates(x), i)


function LineString(points::AbstractVector{LineP{Dim, T, P}}) where {Dim, T, P}
    LineString{Dim, T, P, typeof(points)}(points)
end

"""
    LineString(points::AbstractVector{<: AbstractPoint}, skip = 1)

Creates a LineString from a vector of points
With `skip == 1`, the default, it will connect the line like this:
```julia
points = Point[a, b, c, d]
linestring = LineString(points)
@assert linestring == LineString([a => b, b => c, c => d])
```
"""
function LineString(points::AbstractVector{<: AbstractPoint}, skip = 1)
    LineString(connect(points, LineP, skip))
end

function LineString(points::AbstractVector{<: Pair{P, P}}) where P <: AbstractPoint{N, T} where {N, T}
    LineString(reinterpret(LineP{N, T, P}, points))
end

function LineString(points::AbstractVector{<: AbstractPoint}, faces::AbstractVector{<: LineFace})
    LineString(connect(points, faces))
end

"""
    LineString(points::AbstractVector{<: AbstractPoint}, indices::AbstractVector{<: Integer}, skip = 1)

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
function LineString(points::AbstractVector{<: AbstractPoint}, indices::AbstractVector{<: Integer}, skip = 1)
    faces = connect(indices, LineFace, skip)
    LineString(points, faces)
end

struct Polygon{
        Dim, T <: Real,
        P <: AbstractPoint{Dim, T},
        L <: AbstractVector{<: LineP{Dim, T, P}},
        V <: AbstractVector{L}
    } <: AbstractPolygon{Dim, T}
    exterior::L
    interiors::V
end

Base.:(==)(a::Polygon, b::Polygon) = (a.exterior == b.exterior) && (a.interiors == b.interiors)

function Polygon(exterior::AbstractVector{LineP{Dim, T, P}}, interiors::V) where {Dim, T, P, V <: AbstractVector{<: AbstractVector{LineP{Dim, T, P}}}}
    Polygon{Dim, T, typeof(exterior), V}(exterior, interiors)
end

Polygon(exterior::L) where L <: AbstractVector{<: LineP} = Polygon(exterior, L[])

function Polygon(exterior::AbstractVector{P}, skip::Int = 1) where P <: AbstractPoint{Dim, T} where {Dim, T}
    Polygon(LineString(exterior, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<: Integer}, skip::Int = 1) where P <: AbstractPoint{Dim, T} where {Dim, T}
    Polygon(LineString(exterior, faces, skip))
end

function Polygon(exterior::AbstractVector{P}, faces::AbstractVector{<: LineFace}) where P <: AbstractPoint{Dim, T} where {Dim, T}
    Polygon(LineString(exterior, faces))
end


struct MultiPolygon{
        Dim, T <: Real,
        Element <: AbstractPolygon{Dim, T},
        A <: AbstractVector{Element}
    } <: AbstractVector{Element}

    polygons::A
end

function MultiPolygon(polygons::AbstractVector{P}; meta...) where P <: AbstractPolygon{Dim, T} where {Dim, T}
    nt = values(meta)
    n = keys(nt)
    typs = Tuple{eltype.(values(nt))...}
    m = StructArray(nt)
    PType = MetaPolygon{Dim, T, P, n, typs}
    mpolys = StructArray{PType}((polygons, m))
    MultiPolygon{Dim, T, PType, typeof(mpolys)}(mpolys)
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)

struct MultiLineString{
        Dim, T <: Real,
        Element <: LineString{Dim, T},
        A <: AbstractVector{Element}
    } <: AbstractVector{Element}

    polygons::A
end

struct Mesh{
        Dim, T <: Real,
        Element <: Polytope{Dim, T},
        V <: AbstractVector{Element}
    } <: AbstractVector{Element}

    simplices::V
end

Base.size(x::Mesh) = size(x.simplices)
Base.getindex(x::Mesh, i::Integer) = x.simplices[i]


function Mesh(elements::AbstractVector{<: Polytope{Dim, T}}) where {Dim, T}
    Mesh{Dim, T, eltype(elements), typeof(elements)}(elements)
end

function Mesh(points::AbstractVector{<: AbstractPoint}, faces::AbstractVector{<: AbstractFace})
    Mesh(connect(points, faces))
end

function Mesh(
        points::AbstractVector{<: AbstractPoint}, faces::AbstractVector{<: Integer},
        facetype = TriangleFace, skip = 1
    )
    Mesh(connect(points, connect(faces, facetype, skip)))
end
