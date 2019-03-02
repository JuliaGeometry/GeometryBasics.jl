module GeometryBasics

using StaticArrays, Tables, StructArrays
import GeometryTypes
import Base: convert
"""
Geometry made of N connected points. Connected as one flat geometry, it makes a Line, Ngon / Polygon.
Connected as volume it will be a Simplex / Tri / Cube.
Note That `Polytype{N} where N == 3` denotes a Triangle both as a Simplex or Ngon.
"""
abstract type Polytope{Dimension, T <: Real} end
abstract type AbstractPolygon{Dimension, T} <: Polytope{Dimension, T} end

abstract type AbstractPoint{Dimension, T} <: StaticVector{Dimension, T} end
abstract type AbstractFace{Dimension, T} <: StaticVector{Dimension, T} end
abstract type AbstractSimplex{Dimension, N, T} <: StaticVector{Dimension, T} end

GeometryTypes.@fixed_vector Point AbstractPoint

"""
Fixed Size Polygon, e.g.
N 1-2 : Illegal!
N = 3 : Triangle
N = 4 : Quadrilateral (or Quad, Or tetragon)
N = 5 : Pentagon
...
"""
struct Ngon{
        Dimension, T <: Real,
        N,
        Point <: AbstractPoint{Dimension, T}
    } <: AbstractPolygon{Dimension, T}

    points::SVector{N, Point}
end
# Simplex{D, T, 3} & Ngon{D, T, 3} are both representing a triangle.
# Since Ngon is supposed to be flat and a triangle is flat, lets prefer Ngon
# for triangle:
const Triangle{Dimension, T} = Ngon{Dimension, T, 3, P} where P <: AbstractPoint{Dimension, T}
const Quadrilateral{Dimension, T} = Ngon{Dimension, T, 4, P} where P <: AbstractPoint{Dimension, T}

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
        Dimension, T <: Real,
        N,
        Point <: AbstractPoint{Dimension, T},
    } <: Polytope{Dimension, T}

    points::SVector{N, Point}
end
const Line{Dimension, T, P <: AbstractPoint{Dimension, T}} = Simplex{Dimension, T, 2, P}
const Tetrahedron{Dimension, T, P <: AbstractPoint{Dimension, T}} = Simplex{Dimension, T, 2, P}


struct LineString{
        Dimension, T <: Real,
        V <: AbstractVector{<: AbstractPoint{Dimension, T}}
    } <: AbstractVector{Line{Dimension, T}}

    points::V
end

struct Polygon{
        Dimension, T <: Real,
        L <: AbstractVector{<: Line{Dimension, T}},
        V <: AbstractVector{L}
    } <: AbstractPolygon{Dimension, T}

    exterior::L
    interiors::V
end
function Polygon(exterior::AbstractVector{Line{Dimension, T, P}}, interiors::V) where {Dimension, T, P, V <: AbstractVector{<: AbstractVector{Line{Dimension, T, P}}}}
    Polygon{Dimension, T, typeof(exterior), V}(exterior, interiors)
end

Polygon(exterior::L) where L <: AbstractVector{<: Line} = Polygon(exterior, L[])
function Polygon(exterior::AbstractVector{P}) where {Dim, T, P <: AbstractPoint{Dim, T}}
    Polygon(reinterpret(Line{Dim, T, P}, exterior))
end


struct MetaPolygon{
        Dimension, T <: Real,
        P <: AbstractPolygon{Dimension, T},
        Names, Types
    } <: AbstractPolygon{Dimension, T}

    polygon::P
    meta::NamedTuple{Names, Types}
end

MetaPolygon(polygon::AbstractPolygon; meta...) = MetaPolygon(polygon, values(meta))
function MetaPolygon(polygon::AbstractVector; meta...)
    MetaPolygon(Polygon(polygon); meta...)
end
function MetaPolygon(exterior::L, interior::AbstractVector{L}; meta...) where L
    MetaPolygon(Polygon(exterior, interior); meta...)
end

function Tables.schema(::AbstractVector{MetaPolygon{Dim, T, P, Names, Types}}) where {Dim, T, P, Names, Types}
    Tables.Schema((:polygon, Names...), (P, Types...))
end
function Base.getproperty(x::MetaPolygon{Dim, T, P, Names, Types}, field::Symbol) where {Dim, T, P, Names, Types}
    field === :polygon && return getfield(x, :polygon)
    Base.sym_in(field, Names) && return getfield(getfield(x, :meta), field)
    error("Field $field not part of Element")
end


struct MultiPolygon{
        Dimension, T <: Real,
        Element <: AbstractPolygon{Dimension, T},
        A <: AbstractVector{Element}
    } <: AbstractVector{Element}

    polygons::A
end
Tables.rows(x::MultiPolygon) = x.polygons
Tables.istable(::Type{<:MultiPolygon}) = true
Tables.istable(::Type{<:MultiPolygon{<:Tuple}}) = false
Tables.rowaccess(::Type{<:MultiPolygon}) = true
# Tables.columnaccess(::Type{<:MultiPolygon}) = true

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
        Dimension, T <: Real,
        Element <: LineString{Dimension, T},
        A <: AbstractVector{Element}
    } <: AbstractVector{Element}

    polygons::A
end

struct Mesh{
        Dimension, T <: Real,
        Element <: Polytope{Dimension, T},
        V <: AbstractVector{Element}
    } <: AbstractVector{Element}

    simplices::V
end

"""
FaceView enables to link one array of points via a face array, to generate one
abstract array of elements.
E.g., this becomes possible:
```
x = FaceView(rand(Point3f0, 10), TriangleFace[(1, 2, 3), (2,4, 5), ...])
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


end # module
