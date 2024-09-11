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

"""
    AbstractFace{N, T} <: StaticVector{N, T}

Parent type for all faces. You should inherit from one of the child types instead.
"""
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

struct FaceView{
        T, AVT <: AbstractVector{T}, 
        FVT <: AbstractVector{ <: AbstractFace}
    }

    data::AVT
    faces::FVT
end

const VertexAttributeType{T} = Union{FaceView{T}, AbstractVector{T}}

function Base.vcat(a::FaceView, b::FaceView)
    N = length(a.data)
    return FaceView(
        vcat(a.data, b.data),
        vcat(a.faces, map(f -> typeof(f)(f .+ N), b.faces))
    )
end

"""
    Mesh <: AbstractMesh{Element}
The concrete AbstractMesh type.
"""
struct Mesh{
        Dim, T <: Real, # TODO: Number?
        FT <: AbstractFace,
        Names,
        VAT <: Tuple{AbstractVector{Point{Dim, T}}, Vararg{VertexAttributeType}},
        FVT <: AbstractVector{FT}
    } <: AbstractMesh{Dim, T}

    vertex_attributes::NamedTuple{Names, VAT}
    connectivity::FVT
    views::Vector{UnitRange{Int}}

    function Mesh(
            va::NamedTuple{Names, VAT}, 
            faces::FVT,
            views::Vector{UnitRange{Int}} = UnitRange{Int}[]
        ) where {
            D, T, FT <: AbstractFace, Names,
            VAT <: Tuple{AbstractVector{Point{D, T}}, Vararg{VertexAttributeType}},
            FVT <: AbstractVector{FT}
        }

        # verify type
        if first(Names) !== :position
            error("The first vertex attribute should be a 'position' but is a '$(first(Names))'.")
        end

        if :normals in Names
            @warn "`normals` as a vertex attribute name has been deprecated in favor of `normal` to bring it in line with mesh.position and mesh.uv"
            names = map(name -> ifelse(name === :normals, :normal, name), Names)
            va = NamedTuple{names}(values(va))
        end

        if va.position isa FaceView
            if faces != va.position.faces
                error("position faces do not match gives faces")
            end
            va = NamedTuple(map(pairs(va)) do (key, val)
                return key => (key == :position ? val.data : val)
            end)
        end

        # verify matching faces between face views
        for (name, attrib) in pairs(va)
            name == :position && continue

            if attrib isa FaceView
                if length(attrib.faces) != length(faces)
                    error("Number of faces defined for $name $(length(attrib.faces)) does not match position $(length(faces))")
                end
                
                for (i, (f1, f2)) in enumerate(zip(attrib.faces, faces))
                    if length(f1) != length(f2)
                        error("Length of face $name[$i] = $(length(f1)) does not match position[$i] = $(length(f2))")
                    end
                end
            end
        end

        return new{D, T, FT, keys(va), VAT, FVT}(va, faces, views)
    end
end

# TODO: what shorthands could be useful?
# const GLMesh{D, T} = Mesh{D, T, GLTriangleFace}
# const SimpleGLMesh{D, T} = Mesh{D, T, GLTriangleFace, (:position,)}

@inline function Base.hasproperty(mesh::Mesh, field::Symbol)
    if field === :normals
        @warn "mesh.normals has been deprecated in favor of mesh.normal to bring it in line with mesh.position and mesh.uv"
        return hasfield(mesh, :normal)
    end
    return hasproperty(getfield(mesh, :vertex_attributes), field) || hasfield(Mesh, field)
end
@inline function Base.getproperty(mesh::Mesh, field::Symbol)
    if hasfield(Mesh, field)
        return getfield(mesh, field)
    elseif field === :normals
        @warn "mesh.normals has been deprecated in favor of mesh.normal to bring it in line with mesh.position and mesh.uv"
        return getproperty(mesh, :normal)
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
# function Base.getindex(mesh::Mesh{D, T, <: AbstractVertexFace}, f::AbstractVertexFace) where {D, T}
#     return getfield(mesh, :vertex_attributes).position[f]
# end
# function Base.getindex(::Mesh, f::AbstractMultiFace)
#     error("TODO")
# end

function Base.:(==)(a::Mesh, b::Mesh)
    return (a.vertex_attributes == b.vertex_attributes) && 
           (faces(a) == faces(b)) && (a.views == b.views)
end

function Base.iterate(mesh::Mesh, i=1)
    return i - 1 < length(mesh) ? (mesh[i], i + 1) : nothing
end

# function Mesh(faces::AbstractVector{<:AbstractFace}; views = UnitRange{Int}[], attributes...)
#     return Mesh(NamedTuple(attributes), faces, views)
# end

function Mesh(points::AbstractVector{Point{Dim, T}},
              faces::AbstractVector{<:AbstractFace}; 
              views = UnitRange{Int}[], kwargs...) where {Dim, T}
    return Mesh((position = points, kwargs...), faces, views)
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
Base.keys(mesh::MetaMesh) = keys(getfield(mesh, :meta))

coordinates(mesh::MetaMesh) = coordinates(Mesh(mesh))
faces(mesh::MetaMesh) = faces(Mesh(mesh))
normals(mesh::MetaMesh) = normals(Mesh(mesh))
texturecoordinates(mesh::MetaMesh) = texturecoordinates(Mesh(mesh))
vertex_attributes(mesh::MetaMesh) = vertex_attributes(Mesh(mesh))

meta(@nospecialize(m)) = NamedTuple()
meta(mesh::MetaMesh) = getfield(mesh, :meta)
Mesh(mesh::MetaMesh) = getfield(mesh, :mesh)
Mesh(mesh::Mesh) = mesh


# Shorthand types
# used in meshes.jl
const SimpleMesh{N, T, FT} = Mesh{N, T, FT, (:position,), Tuple{Vector{Point{N, T}}}, Vector{FT}}
const SimpleTriangleMesh{N} = SimpleMesh{N, Float32, GLTriangleFace}
#
const NormalMesh{N, T, FT}   = Mesh{N, T, FT, (:position, :normal), Tuple{Vector{Point{N, T}}, Vector{Vec3f}}, Vector{FT}}
const NormalUVMesh{N, T, FT} = Mesh{N, T, FT, (:position, :normal, :uv), Tuple{Vector{Point{N, T}}, Vector{Vec3f}, Vector{Vec2f}}, Vector{FT}}
const GLNormalMesh{N}   = NormalMesh{N, Float32, GLTriangleFace}
const GLNormalUVMesh{N} = NormalUVMesh{N, Float32, GLTriangleFace}