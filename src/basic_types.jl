"""
Abstract Geometry in R{Dim} with Number type T
"""
abstract type AbstractGeometry{Dim,T<:Number} end
abstract type GeometryPrimitive{Dim,T} <: AbstractGeometry{Dim,T} end
Base.ndims(::AbstractGeometry{Dim}) where {Dim} = Dim

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

@propagate_inbounds function Base.getindex(points::AbstractVector{P}, face::F) where {P<: Point, F <: AbstractFace}
    Polytope(P, F)(map(i-> points[i], face.data))
end

@fixed_vector SimplexFace AbstractSimplexFace

const TetrahedronFace{T} = SimplexFace{4,T}
Face(::Type{<:SimplexFace{N}}, ::Type{T}) where {N,T} = SimplexFace{N,T}

@fixed_vector NgonFace AbstractNgonFace

const LineFace{T} = NgonFace{2,T}
const TriangleFace{T} = NgonFace{3,T}
const QuadFace{T} = NgonFace{4,T}
const GLTriangleFace = TriangleFace{GLIndex}

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
function Polytope(::Type{Point{Dim,T}},
                  ::Type{<:AbstractNgonFace{N,IT}}) where {N,Dim,T,IT}
    return Ngon{Dim,T,N}
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
const GLTriangleElement = Triangle{3,Float32}

faces(x::Ngon{Dim, T, N}) where {Dim, T, N} = [NgonFace{N, Int}(ntuple(identity, N))]

Base.show(io::IO, x::Triangle) = print(io, "Triangle(", join(x, ", "), ")")

const Quadrilateral{Dim,T} = Ngon{Dim,T,4}

Base.show(io::IO, x::Quadrilateral) = print(io, "Quadrilateral(", join(x, ", "), ")")

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
function Polytope(::Type{Point{Dim,T}}, ::Type{<:AbstractSimplexFace{N}}) where {N,Dim,T}
    return Simplex{Dim,T,N}
end

"""
The fully concrete Simplex type, when constructed from a point type!
"""
function Polytope(::Type{<:NSimplex{N}}, P::Type{Point{NDim,T}}) where {N,NDim,T}
    return Simplex{NDim,T,N}
end
Base.show(io::IO, x::Line) = print(io, "Line(", x[1], " => ", x[2], ")")


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
    # if we just map over it, it won't infer the element type correctly!
    int = typeof(exterior)[]
    foreach(x -> push!(int, x), interior)
    return Polygon(exterior, int)
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
    AbstractMesh

An abstract mesh is a collection of Polytope elements (Simplices / Ngons).
The connections are defined via faces(mesh), the coordinates of the elements are returned by
coordinates(mesh). Arbitrary meta information can be attached per point or per face
"""
abstract type AbstractMesh{Dim, T} <: AbstractGeometry{Dim, T} end

"""
    Mesh <: AbstractMesh{Element}
The concrete AbstractMesh type.
"""
struct Mesh{Dim, T<:Number, V<:AbstractVector{Point{Dim, T}}, C <: AbstractVector{<: AbstractFace}} <: AbstractMesh{Dim, T}
    vertices::V
    connectivity::C
end
coordinates(mesh::Mesh) = mesh.vertices
faces(mesh::Mesh) = mesh.connectivity
Base.getindex(mesh::Mesh, i::Integer) = mesh.vertices[mesh.connectivity[i]]
Base.length(mesh::Mesh) = length(mesh.connectivity)
Base.:(==)(a::Mesh, b::Mesh) = coordinates(a) == coordinates(b) && faces(a) == faces(b)

function Base.iterate(mesh::Mesh, i=1)
    return i - 1 < length(mesh) ? (mesh[i], i + 1) : nothing
end

function Mesh(points::AbstractVector{Point{Dim, T}},
              faces::AbstractVector{<:AbstractFace}) where {Dim, T}
    return Mesh{Dim, T, }(points, faces)
end

function Mesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(points, connect(faces, facetype, skip))
end

struct MetaMesh{Dim, T, M <: AbstractMesh{Dim, T}, Names, Types} <: AbstractMesh{Dim, T}
    mesh::M
    meta::NamedTuple{Names, Types}
    function MetaMesh(mesh::AbstractMesh{Dim, T}, meta::NamedTuple{Names, Types}) where {Dim, T, Names, Types}
        new{Dim, T, typeof(mesh), Names, Types}(mesh, meta)
    end
end

function MetaMesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractFace}; meta...)
    MetaMesh(Mesh(points, faces), values(meta))
end

function MetaMesh(m::AbstractMesh; meta...)
    MetaMesh(m, values(meta))
end

@inline Base.hasproperty(mesh::MetaMesh, field::Symbol) = hasproperty(getfield(mesh, :meta), field)
@inline Base.getproperty(mesh::MetaMesh, field::Symbol) = getproperty(getfield(mesh, :meta), field)
@inline Base.propertynames(mesh::MetaMesh) = propertynames(getfield(mesh, :meta))

coordinates(mesh::MetaMesh) = coordinates(Mesh(mesh))
faces(mesh::MetaMesh) = faces(Mesh(mesh))
normals(mesh::MetaMesh) = hasproperty(mesh, :normals) ? mesh.normals : nothing
texturecoordinates(mesh::MetaMesh) = hasproperty(mesh, :uv) ? mesh.uv : nothing

meta(mesh::MetaMesh) = getfield(mesh, :meta)
Mesh(mesh::MetaMesh) = getfield(mesh, :mesh)
