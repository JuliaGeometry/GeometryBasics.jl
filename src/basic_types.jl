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

abstract type AbstractSimplex{Dim,T} <: Polytope{Dim,T} end

@propagate_inbounds function Base.getindex(points::AbstractVector{P}, face::F) where {P<: Point, F <: AbstractFace}
    return Polytope(P, F)(map(i-> points[i], face.data))
end

@propagate_inbounds function Base.getindex(elements::AbstractVector, face::F) where {F <: AbstractFace}
    return map(i-> elements[i], face.data)
end

@fixed_vector SimplexFace = AbstractSimplexFace

const TetrahedronFace{T} = SimplexFace{4,T}
Face(::Type{<:SimplexFace{N}}, ::Type{T}) where {N,T} = SimplexFace{N,T}

"""
Face index, connecting points to form an Ngon
"""

@fixed_vector NgonFace = AbstractNgonFace

const LineFace{T} = NgonFace{2,T}
const TriangleFace{T} = NgonFace{3,T}
const QuadFace{T} = NgonFace{4,T}
const GLTriangleFace = TriangleFace{GLIndex}

function Base.show(io::IO, x::TriangleFace{T}) where {T}
    return print(io, "TriangleFace(", join(x, ", "), ")")
end

Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F

struct MultiFace{N, T, FaceType <: AbstractFace{N, T}, Names, M} <: AbstractFace{N, T}
    faces::NamedTuple{Names, NTuple{M, FaceType}}

    function MultiFace(nt::NamedTuple{Names, NTuple{M, FT}}) where {N, T, FT <: AbstractFace{N, T}, Names, M}
        if FT <: MultiFace
            error("A MultiFace cannot contain MultiFaces.")
        end

        return new{N, T, FT, Names, M}(nt)
    end
end

MultiFace(; kwargs...) = MultiFace(NamedTuple(kwargs))
MultiFace{Names}(args...) where {Names} = MultiFace(NamedTuple{Names}(args))
MultiFace{Names}(args::Tuple{Vararg{<: AbstractFace}}) where {Names} = MultiFace(NamedTuple{Names}(args))
MultiFace{Names, FT}(args) where {Names, FT <: AbstractFace} = MultiFace(NamedTuple{Names}(FT.(args)))

function MultiFace{Names}(f::MultiFace) where {Names}
    return MultiFace{Names}(getproperty.((f,), Names))
end

Base.getindex(f::MultiFace, i::Integer) = Base.getindex(getfield(f, :faces), i)
@inline Base.hasproperty(f::MultiFace, field::Symbol) = hasproperty(getfield(f, :faces), field)
@inline Base.getproperty(f::MultiFace, field::Symbol) = getproperty(getfield(f, :faces), field)
@inline Base.propertynames(f::MultiFace) = propertynames(getfield(f, :faces))
@inline Base.propertynames(::Type{<: MultiFace{N, T, FT, Names}}) where {N, T, FT, Names} = Names
Base.eltype(::MultiFace{N, T, FT}) where {N, T, FT} = FT
Base.eltype(::Type{<: MultiFace{N, T, FT}}) where {N, T, FT} = FT

function simplify_faces(::Type{MF1}, fs::AbstractVector{MF2}) where {MF1 <: MultiFace, MF2 <: MultiFace}
    return simplify_faces(propertynames(MF1), fs)
end

function simplify_faces(names::NTuple{N, Symbol}, fs::AbstractVector{MF2}) where {N, MF2 <: MultiFace}
    return map(f -> MultiFace{names}(f), fs)
end

# TODO: Some shorthands
const NormalFace = MultiFace{(:position, :normal)}
const NormalUVFace = MultiFace{(:position, :normal, :uv)}

# TODO: enable something like NormalUVFace{QuadFace}[...]

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
struct Simplex{Dim,T<:Real,N} <: AbstractSimplex{Dim,T}
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
struct Polygon{Dim,T<:Real} <: AbstractPolygon{Dim,T}
    exterior::Vector{Point{Dim, T}}
    interiors::Vector{Vector{Point{Dim, T}}}
end

Base.copy(x::Polygon) = Polygon(copy(x.exterior), deepcopy(x.interiors))

function Base.:(==)(a::Polygon, b::Polygon)
    return (a.exterior == b.exterior) && (a.interiors == b.interiors)
end

function Polygon(
        exterior::AbstractVector{Point{Dim,T}},
        interiors::AbstractVector{AbstractVector{Point{Dim,T}}}) where {Dim, T}
    tov(v) = convert(Vector{Point{Dim, T}}, v)
    return Polygon{Dim,T}(tov(exterior), map(tov, interiors))
end

Polygon(exterior::AbstractVector{Point{N, T}}) where {N, T} = Polygon(exterior, Vector{Point{N, T}}[])

function Polygon(exterior::AbstractVector{Point{Dim,T}},
                 faces::AbstractVector{<:LineFace}) where {Dim, T}
    return Polygon(LineString(exterior, faces))
end

function Polygon(exterior::AbstractGeometry{Dim, T}, interior::AbstractVector=[]) where {Dim, T}
    to_p(v) = decompose(Point{Dim, T}, v)
    int = Vector{Point{Dim, T}}[]
    for i in interior
        push!(int, to_p(i))
    end
    return Polygon(to_p(exterior), int)
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
struct MultiPolygon{Dim, T<:Real} <: AbstractGeometry{Dim, T}
    polygons::Vector{<:AbstractPolygon{Dim,T}}
end

function MultiPolygon(polygons::AbstractVector{<:AbstractPolygon{Dim,T}}) where {Dim,T}
    return MultiPolygon(convert(Vector{eltype(polygons)}, polygons))
end

Base.getindex(mp::MultiPolygon, i) = mp.polygons[i]
Base.size(mp::MultiPolygon) = size(mp.polygons)
Base.length(mp::MultiPolygon) = length(mp.polygons)

"""
    LineString(points::AbstractVector{<:Point})
A LineString is a geometry of connected line segments
"""
struct LineString{Dim, T<:Real} <: AbstractGeometry{Dim, T}
    points::Vector{Point{Dim, T}}
end
Base.length(ls::LineString) = length(coordinates(ls))
coordinates(ls::LineString) = ls.points

struct MultiLineString{Dim, T<:Real} <: AbstractGeometry{Dim, T}
    linestrings::Vector{LineString{Dim, T}}
end

function MultiLineString(linestrings::AbstractVector{L}) where {L<:LineString}
    return MultiLineString(convert(Vector{L}, linestrings))
end

Base.getindex(ms::MultiLineString, i) = ms.linestrings[i]
Base.size(ms::MultiLineString) = size(ms.linestrings)
Base.length(mpt::MultiLineString) = length(mpt.linestrings)

"""
    MultiPoint(points::AbstractVector{AbstractPoint})

A collection of points
"""
struct MultiPoint{Dim,T<:Real} <: AbstractGeometry{Dim, T}
    points::Vector{Point{Dim, T}}
end

function MultiPoint(points::AbstractVector{Point{Dim, T}}) where {Dim,T}
    return MultiPoint(convert(Vector{Point{Dim, T}}, points))
end

Base.getindex(mpt::MultiPoint, i) = mpt.points[i]
Base.size(mpt::MultiPoint) = size(mpt.points)
Base.length(mpt::MultiPoint) = length(mpt.points)

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
struct Mesh{
        Dim, T <: Real, # TODO: Number?
        FaceType <: AbstractFace,
        Names,
        VertexAttribTypes <: Tuple{AbstractVector{Point{Dim, T}}, Vararg{AbstractVector}},
        FaceVecType <: AbstractVector{FaceType}
    } <: AbstractMesh{Dim, T}

    vertex_attributes::NamedTuple{Names, VertexAttribTypes}
    connectivity::FaceVecType
    views::Vector{UnitRange{Int}}

    function Mesh(
            va::NamedTuple{Names, VAT}, 
            f::FVT, 
            views::Vector{UnitRange{Int}} = UnitRange{Int}[]
        ) where {
            D, T, FT, Names,
            VAT <: Tuple{AbstractVector{Point{D, T}}, Vararg{AbstractVector}},
            FVT <: AbstractVector{FT}
        }

        # verify type
        if first(Names) !== :position
            error("The first vertex attribute should be a 'position' but is a '$(first(Names))'.")
        end

        if FT <: MultiFace
            f_names = propertynames(FT)
            # if Names != f_names
            #     error(
            #         "Cannot construct a mesh with vertex attribute names $Names and MultiFace " * 
            #         "attribute names $f_names. These must include the same names in the same order."
            #     )
            # end
            if Names == f_names
                # all good
            elseif issubset(Names, f_names)
                # remove the extras/fix order
                f = simplify_faces(Names, f)
            else
                error(
                    "Cannot construct a mesh with vertex attribute names $Names and MultiFace " * 
                    "attribute names $f_names. These must include the same names in the same order."
                )
            end
        elseif MultiFace <: FT 
            # TODO: This is supposed to catch mixed types like
            # [MultiFace(position = f1, normal = f2), MultiFace(position = f3)]
            # but really just catches AbstractFace{N, T}[]. Technically we can
            # probably handle mixtures of MultiFace and other Face types, but do
            # we want to bother? Also do we want to allow mixtures of e.g.
            # TriangleFace and QuadFace?
            error("Face vectors that may include `MultiFace`s with different names are not allowed. (Type $FT too abstract.)")
        end

        return new{D, T, eltype(typeof(f)), Names, VAT, typeof(f)}(va, f, views)
    end
end

@inline function Base.hasproperty(mesh::Mesh, field::Symbol)
    return hasproperty(getfield(mesh, :vertex_attributes), field) || hasfield(Mesh, field)
end
@inline function Base.getproperty(mesh::Mesh, field::Symbol)
    if hasfield(Mesh, field)
        return getfield(mesh, field)
    else
        return getproperty(getfield(mesh, :vertex_attributes), field)
    end
end
@inline function Base.propertynames(mesh::Mesh)
    return (fieldnames(Mesh)..., propertynames(getfield(mesh, :vertex_attributes))...)
end

coordinates(mesh::Mesh) = mesh.position
faces(mesh::Mesh) = mesh.connectivity
normals(mesh::Mesh) = hasproperty(mesh, :normal) ? mesh.normal : nothing
texturecoordinates(mesh::Mesh) = hasproperty(mesh, :uv) ? mesh.uv : nothing
vertex_attributes(mesh::Mesh) = getfield(mesh, :vertex_attributes)

Base.getindex(mesh::Mesh, i::Integer) = mesh[mesh.connectivity[i]]
Base.length(mesh::Mesh) = length(mesh.connectivity)

# TODO: temp
function Base.getindex(mesh::Mesh{D, T, <: AbstractFace, (:position,)}, f::AbstractFace) where {D, T}
    return getfield(mesh, :vertex_attributes).position[f]
end
function Base.getindex(::Mesh, f::MultiFace)
    error("TODO")
end
function Base.getindex(::Mesh, f::AbstractFace)
    error("TODO")
end

function Base.:(==)(a::Mesh, b::Mesh)
    return a.vertex_attributes == b.vertex_attributes && 
           faces(a) == faces(b) && a.views == b.views
end

function Base.iterate(mesh::Mesh, i=1)
    return i - 1 < length(mesh) ? (mesh[i], i + 1) : nothing
end

function Mesh(faces::AbstractVector{<:AbstractFace}; attributes...)
    return Mesh(NamedTuple(attributes), faces)
end

function Mesh(points::AbstractVector{Point{Dim, T}},
              faces::AbstractVector{<:AbstractFace}; kwargs...) where {Dim, T}
    return Mesh((position = points, kwargs...), faces)
end

function Mesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(points, connect(faces, facetype, skip))
end


struct MetaMesh{Dim, T, M <: AbstractMesh{Dim, T}} <: AbstractMesh{Dim, T}
    mesh::M
    meta::Dict{Symbol, Any}
end

function MetaMesh(mesh::AbstractMesh; kwargs...)
    MetaMesh(mesh, Dict{Symbol, Any}(kwargs))
end

function MetaMesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractFace}; kwargs...)
    MetaMesh(Mesh(points, faces), Dict{Symbol, Any}(kwargs))
end

# TODO: Do we want to access meta here?
@inline function Base.hasproperty(mesh::MetaMesh, field::Symbol)
    return hasfield(MetaMesh, field) || hasproperty(getfield(mesh, :mesh), field)
end
@inline function Base.getproperty(mesh::MetaMesh, field::Symbol)
    if hasfield(MetaMesh, field)
        return getfield(mesh, field)
    else
        return getproperty(getfield(mesh, :mesh), field)
    end
end
@inline function Base.propertynames(mesh::MetaMesh)
    return (fieldnames(MetaMesh)..., propertynames(getfield(mesh, :mesh))...)
end

# TODO: or via getindex?
Base.haskey(mesh::MetaMesh, key::Symbol) = haskey(getfield(mesh, :meta), key)
Base.get(f, mesh::MetaMesh, key::Symbol) = get(f, getfield(mesh, :meta), key)
Base.get!(f, mesh::MetaMesh, key::Symbol) = get!(f, getfield(mesh, :meta), key)
Base.getindex(mesh::MetaMesh, key::Symbol) = getindex(getfield(mesh, :meta), key)
Base.setindex!(mesh::MetaMesh, value, key::Symbol) = setindex!(getfield(mesh, :meta), value, key)
Base.delete!(mesh::MetaMesh, key::Symbol) = delete!(getfield(mesh, :meta), key)

coordinates(mesh::MetaMesh) = coordinates(Mesh(mesh))
faces(mesh::MetaMesh) = faces(Mesh(mesh))
normals(mesh::MetaMesh) = normals(Mesh(mesh))
texturecoordinates(mesh::MetaMesh) = texturecoordinates(Mesh(mesh))
vertex_attributes(mesh::MetaMesh) = vertex_attributes(Mesh(mesh))

meta(@nospecialize(m)) = NamedTuple()
meta(mesh::MetaMesh) = getfield(mesh, :meta)
Mesh(mesh::MetaMesh) = getfield(mesh, :mesh)
Mesh(mesh::Mesh) = mesh
