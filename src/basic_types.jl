"""
    abstract type AbstractGeometry{Dimension, T<:Number}

Base type for geometry types like GeometryPrimites and Polytopes.
"""
abstract type AbstractGeometry{Dim,T<:Number} end
abstract type GeometryPrimitive{Dim,T} <: AbstractGeometry{Dim,T} end
Base.ndims(::AbstractGeometry{Dim}) where {Dim} = Dim

"""
    Polytope{Dim, T} <: AbstractGeometry{Dim, T}

A Polytope is the generalization of a Polygon to higher dimensions, i.e. a 
geometric object consisting of flat faces.

A `Polygon` and `Ngon` are both 2D `Polytope`s. A `Simplex` is the simplest
`Polytope` in arbitrary dimensions.
"""
abstract type Polytope{Dim,T} <: AbstractGeometry{Dim,T} end
abstract type AbstractPolygon{Dim,T} <: Polytope{Dim,T} end

"""
    AbstractFace{N_indices, T} <: StaticVector{N_indices, T}

Parent type for all face types. The standard face type is typically a 
`GLTriangleFace = NgonFace{3, GLIndex}`.
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


@fixed_vector NgonFace = AbstractNgonFace

"""
    NgonFace{N, T}

A planar face connecting N vertices. Shorthands include:
- `LineFace{T} = NgonFace{2,T}`
- `TriangleFace{T} = NgonFace{3,T}`
- `QuadFace{T} = NgonFace{4,T}`
- `GLTriangleFace = TriangleFace{GLIndex}`
"""
NgonFace

const LineFace{T} = NgonFace{2,T}
const TriangleFace{T} = NgonFace{3,T}
const QuadFace{T} = NgonFace{4,T}
const GLTriangleFace = TriangleFace{GLIndex}

function Base.show(io::IO, x::NgonFace{N, T}) where {N, T}
    if N == 2
        name = "LineFace{$T}"
    elseif N == 3
        if T == GLIndex
            name = "GLTriangleFace"
        else
            name = "TriangleFace{$T}"
        end
    elseif N == 4
        name = "QuadFace{$T}"
    else
        name = "NgonFace{$N, $T}"
    end

    return print(io, name, "(", join(value.(x), ", "), ")")
end

Face(::Type{<:NgonFace{N}}, ::Type{T}) where {N,T} = NgonFace{N,T}
Face(F::Type{NgonFace{N,FT}}, ::Type{T}) where {FT,N,T} = F

@propagate_inbounds Base.getindex(x::Polytope, i::Integer) = coordinates(x)[i]
@propagate_inbounds Base.iterate(x::Polytope) = iterate(coordinates(x))
@propagate_inbounds Base.iterate(x::Polytope, i) = iterate(coordinates(x), i)

"""
    Ngon{D, T, N}(points::NTuple{N, Point{D, T}})

Defines a flat polygon (without holes) in D dimensional space using N points, e.g.:
- N 1-2 : Illegal!
- N = 3 : Triangle
- N = 4 : Quadrilateral (or Quad, Or tetragon)
- N = 5 : Pentagon
- ...

For polygons with holes, see `Polygon`.
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
    Polytope(::Type{<: Point}, ::Type{<: AbstractNgonFace})

The Ngon Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(::Type{Point{Dim,T}},
                  ::Type{<:AbstractNgonFace{N,IT}}) where {N,Dim,T,IT}
    return Ngon{Dim,T,N}
end

"""
    Polytope(::Type{<: Ngon}, P::Type{<: Point})

The fully concrete Ngon type, when constructed from a point type!
"""
function Polytope(::Type{<:Ngon{_D, _T, N}}, P::Type{Point{NDim,T}}) where {N,NDim,T, _D,_T}
    return Ngon{NDim,T,N}
end
function Polytope(::Type{<:Ngon{_D, _T, N} where {_D,_T}}, P::Type{Point{NDim,T}}) where {N,NDim,T}
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
    Simplex{D, T<:Real, N}(points::NTuple{N, Point{D, T}})

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
    Polytope(::Type{Point{Dim,T}}, ::Type{<:AbstractSimplexFace{N}})

The Simplex Polytope element type when indexing an array of points with a SimplexFace
"""
function Polytope(::Type{Point{Dim,T}}, ::Type{<:AbstractSimplexFace{N}}) where {N,Dim,T}
    return Simplex{Dim,T,N}
end

"""
    Polytope(::Type{<:NSimplex{N}}, P::Type{Point{NDim,T}})

The fully concrete Simplex type, when constructed from a point type!
"""
function Polytope(::Type{<:NSimplex{N}}, P::Type{Point{NDim,T}}) where {N,NDim,T}
    return Simplex{NDim,T,N}
end
Base.show(io::IO, x::Line) = print(io, "Line(", x[1], " => ", x[2], ")")


"""
    Polygon(exterior::AbstractVector{<:Point})
    Polygon(exterior::AbstractVector{<:Point}, interiors::Vector{<:AbstractVector{<:Point}})

Constructs a polygon from a set of exterior points. If interiors are given, each
of them cuts away from the Polygon.
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

A collection of polygons
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

A LineString is a collection of points connected by line segments.
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
    FaceView(data, faces)

A FaceView is alternative to passing a vertex attribute directly to a mesh. It 
bundles `data` with a different set of `faces` which replace the default faces
of a mesh. Doing this allows you to have the `data` in a different order from 
vertices in the mesh, potentially avoiding duplication.

You can get the data of a FaceView with `values(faceview)` and the faces with
`faces(faceview)`.
"""
struct FaceView{T, AVT <: AbstractVector{T}, FVT <: AbstractVector{<: AbstractFace}}
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

faces(x::FaceView) = x.faces
Base.values(x::FaceView) = x.data
facetype(x::FaceView) = eltype(x.faces)
Base.getindex(x::FaceView, f::AbstractFace) = getindex(values(x), f)
Base.isempty(x::FaceView) = isempty(values(x))
Base.:(==)(a::FaceView, b::FaceView) = (values(a) == values(b)) && (faces(a) == faces(b))

# TODO: maybe underscore this as it requires care to make sure all FaceViews and
#       mesh faces stay in sync
convert_facetype(::Type{FT}, x::AbstractVector) where {FT <: AbstractFace} = x
function convert_facetype(::Type{FT}, x::FaceView) where {FT <: AbstractFace}
    if eltype(faces(x)) != FT
        return FaceView(values(x), decompose(FT, faces(x)))
    end
    return x
end

function verify(fs::AbstractVector{FT}, fv::FaceView, name = nothing) where {FT <: AbstractFace}
    if isconcretetype(FT) && (FT == facetype(fv))
        return true
    end

    if length(faces(fv)) != length(fs)
        error("Number of faces given in FaceView $(length(faces(fv))) does not match reference $(length(fs))")
    end
    
    for (i, (f1, f2)) in enumerate(zip(faces(fv), fs))
        if length(f1) != length(f2)
            error("Length of face $i = $(length(f1)) does not match reference with $(length(f2))")
        end
    end

    return true
end

# Dodgy definitions... (since attributes can be FaceView or Array it's often 
#                       useful to treat a FaceView like the vertex data it contains)
Base.length(x::FaceView) = length(values(x))
# Base.iterate(x::FaceView) = iterate(values(x))
# Base.getindex(x::FaceView, i::Integer) = getindex(values(x), i)
# Taken from Base/arrayshow.jl
function Base.show(io::IO, ::MIME"text/plain", X::FaceView)
    summary(io, X)
    isempty(X) && return
    print(io, ":")

    if get(io, :limit, false)::Bool && displaysize(io)[1]-4 <= 0
        return print(io, " …")
    else
        println(io)
    end

    io = IOContext(io, :typeinfo => eltype(values(X)))

    recur_io = IOContext(io, :SHOWN_SET => values(X))
    Base.print_array(recur_io, values(X))
end



"""
    AbstractMesh

An abstract mesh is a collection of Polytope elements (Simplices / Ngons). The 
connections are defined via faces(mesh) and the coordinates of the elements are 
returned by coordinates(mesh).
"""
abstract type AbstractMesh{Dim, T} <: AbstractGeometry{Dim, T} end

"""
    Mesh{D, T, FaceType, FaceArrayType} <: AbstractMesh{D, T} <: AbstractGeometry{D, T}

The type of a concrete mesh. It can hold arbitrary vertex data (including FaceView's).
"""
struct Mesh{
        Dim, T <: Real, # TODO: Number?
        FT <: AbstractFace,
        FVT <: AbstractVector{FT}
    } <: AbstractMesh{Dim, T}

    vertex_attributes::Dict{Symbol, VertexAttributeType} 
    faces::FVT
    views::Vector{UnitRange{Int}}

    function Mesh(
            va::Dict{Symbol, VertexAttributeType}, 
            fs::FVT,
            views::Vector{UnitRange{Int}} = UnitRange{Int}[]
        ) where {FT <: AbstractFace, FVT <: AbstractVector{FT}}

        # verify type
        if !haskey(va, :position )
            error("Vertex attributes must have a :position attribute.")
        end

        if haskey(va, :normals)
            @warn "`normals` as a vertex attribute name has been deprecated in favor of `normal` to bring it in line with mesh.position and mesh.uv"
            va[:normal] = pop!(va, :normal)
        end

        # verify that faces of FaceViews match `fs` (in length per face)
        N = maximum(f -> value(maximum(f)), fs, init = 0)
        for (name, attrib) in va
            name == :position && continue

            if attrib isa FaceView
                try
                    verify(fs, attrib)
                catch e
                    @error "Failed to verify $name attribute:"
                    rethrow(e)
                end
            else
                length(attrib) < N && error("Failed to verify $name attribute:\nFaces address $N vertex attributes but only $(length(attrib)) are present.")
            end
        end

        D = length(eltype(va[:position]))
        T = eltype(eltype(va[:position]))

        return new{D, T, FT, FVT}(va, fs, views)
    end
end

@inline function Base.hasproperty(mesh::Mesh, field::Symbol)
    if field === :normals
        @warn "mesh.normals has been deprecated in favor of mesh.normal to bring it in line with mesh.position and mesh.uv"
        return hasproperty(mesh, :normal)
    end
    return haskey(getfield(mesh, :vertex_attributes), field) || hasfield(Mesh, field)
end
@inline function Base.getproperty(mesh::Mesh, field::Symbol)
    if hasfield(Mesh, field)
        return getfield(mesh, field)
    elseif field === :normals
        @warn "mesh.normals has been deprecated in favor of mesh.normal to bring it in line with mesh.position and mesh.uv"
        return getproperty(mesh, :normal)
    else
        return getindex(getfield(mesh, :vertex_attributes), field)
    end
end
@inline function Base.propertynames(mesh::Mesh)
    return (fieldnames(Mesh)..., keys(getfield(mesh, :vertex_attributes))...)
end

coordinates(mesh::Mesh) = mesh.position
faces(mesh::Mesh) = mesh.faces
normals(mesh::Mesh) = hasproperty(mesh, :normal) ? mesh.normal : nothing
texturecoordinates(mesh::Mesh) = hasproperty(mesh, :uv) ? mesh.uv : nothing

"""
    vertex_attributes(mesh::Mesh)

Returns a dictionairy containing the vertex attributes of the given mesh. 
Mutating these will change the mesh.
"""
vertex_attributes(mesh::Mesh) = getfield(mesh, :vertex_attributes)

Base.getindex(mesh::Mesh, i::Integer) = mesh.position[mesh.faces[i]]
Base.length(mesh::Mesh) = length(mesh.faces)

function Base.:(==)(a::Mesh, b::Mesh)
    return (a.vertex_attributes == b.vertex_attributes) && 
           (faces(a) == faces(b)) && (a.views == b.views)
end

function Base.iterate(mesh::Mesh, i=1)
    return i - 1 < length(mesh) ? (mesh[i], i + 1) : nothing
end

function add_vertex_attribute!(mesh::Mesh, val::AbstractVector, name::Symbol)
    @boundscheck begin
        N = maximum(f -> value(maximum(f)), faces(mesh))
        length(val) < N && error("Given vertex data not large enough to be adressed by all faces ($N required, $(length(val)) given)")
    end
    mesh.vertex_attributes[name] = val
    return mesh
end

function add_vertex_attribute!(mesh::Mesh, val::FaceView, name::Symbol)
    @boundscheck verify(faces(mesh), val)
    mesh.vertex_attributes[name] = val
    return mesh
end

"""
    Mesh(faces[; views, attributes...])
    Mesh(positions, faces[; views])
    Mesh(positions, faces::AbstractVector{<: Integer}[; facetype = TriangleFace, skip = 1])
    Mesh(; attributes...)

Constructs a mesh from the given arguments.

If `positions` are given explicitly, they are merged with other vertex attributes 
under the name `position`. Otherwise they must be part of `attributes`. If `faces`
are not given `attributes.position` must be a FaceView.

Any other vertex attribute can be either an `AbstractVector` or a `FaceView` 
thereof. Every vertex attribute that is an `AbstractVector` must be sufficiently
large to be indexable by `mesh.faces`. Every vertex attribute that is a `FaceView`
must contain similar faces to `mesh.faces`, i.e. contain the same number of faces
and have faces of matching length.

`views` can be defined optionally to implicitly split the mesh into multi 
sub-meshes. This is done by providing ranges for indexing faces which correspond
to the sub-meshes. By default this is left empty.
"""
function Mesh(faces::AbstractVector{<:AbstractFace}; views::Vector{UnitRange{Int}} = UnitRange{Int}[], attributes...)
    return Mesh(Dict{Symbol, VertexAttributeType}(attributes), faces, views)
end

function Mesh(points::AbstractVector{Point{Dim, T}},
              faces::AbstractVector{<:AbstractFace}; 
              views = UnitRange{Int}[], kwargs...) where {Dim, T}
    va = Dict{Symbol, VertexAttributeType}(kwargs)
    va[:position] = points
    return Mesh(va, faces, views)
end

function Mesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:Integer},
              facetype=TriangleFace, skip=1)
    return Mesh(points, connect(faces, facetype, skip))
end

function Mesh(; kwargs...)
    va = Dict{Symbol, VertexAttributeType}(kwargs)
    fs = faces(va[:position]::FaceView)
    va[:position] = values(va[:position])
    return Mesh(va, fs)
end

struct MetaMesh{Dim, T, M <: AbstractMesh{Dim, T}} <: AbstractMesh{Dim, T}
    mesh::M
    meta::Dict{Symbol, Any}
end


"""
    MetaMesh(mesh; metadata...)
    MetaMesh(positions, faces; metadata...)

Constructs a MetaMesh either from another `mesh` or by constructing another mesh
with the given `positions` and `faces`. Any keyword arguments given will be
stored in the `meta` field in `MetaMesh`.

This struct is meant to be used for storage of non-vertex data. Any vertex 
related data should be stored as a vertex attribute in `Mesh`. 

The metadata added to the MetaMesh can be manipulated with Dict-like operations 
(getindex, setindex!, get, delete, keys, etc). Vertex attributes can be accessed 
via fields and the same getters as mesh. The mesh itself can be retrieved with 
`Mesh(metamesh)`.
"""
function MetaMesh(mesh::AbstractMesh; kwargs...)
    MetaMesh(mesh, Dict{Symbol, Any}(kwargs))
end

function MetaMesh(points::AbstractVector{<:Point}, faces::AbstractVector{<:AbstractFace}; kwargs...)
    MetaMesh(Mesh(points, faces), Dict{Symbol, Any}(kwargs))
end


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