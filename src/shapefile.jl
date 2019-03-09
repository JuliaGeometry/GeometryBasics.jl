import Base.==
import GeoInterface
using GeometryTypes


struct ShapeRect{T}
    left::T
    bottom::T
    right::T
    top::T
end

struct NullShape <: GeoInterface.AbstractGeometry
end

struct Interval{T}
    left::T
    right::T
end

struct MultiLineString {N, T <: Real, V <: AbstractVector{<: AbstractPoint{N, T}}} <: AbstractVector{Line{N, T}}
    points::V
end


struct Polyline{N, T, A <: AbstractVector{<: AbstractVector{ <: AbstractPoint{N, T}}}}
    points::A
end

struct MeasureMeta{T, M}
    data::T
    measures::M
end

struct Polygon{T} <: GeoInterface.AbstractMultiPolygon
    parts::Vector{Int32}
    points::Vector{Point{T}}
end

Base.show(io::IO,p::Polygon{T}) where {T} = print(io,"Polygon(",length(p.points)," ",T," Points)")

struct PolygonM{T,M} <: GeoInterface.AbstractMultiPolygon
    MBR::Rect{T}
    parts::Vector{Int32}
    points::Vector{Point{T}}
    measures::Vector{M}
end

struct PolygonZ{T,M} <: GeoInterface.AbstractMultiPolygon
    MBR::Rect{T}
    parts::Vector{Int32}
    points::Vector{Point{T}}
    zvalues::Vector{T}
    measures::Vector{M}
end

struct MultiPoint{T} <: GeoInterface.AbstractMultiPoint
    MBR::Rect{T}
    points::Vector{Point{T}}
end

struct MultiPointM{T,M} <: GeoInterface.AbstractMultiPoint
    MBR::Rect{T}
    points::Vector{Point{T}}
    measures::Vector{M}
end

struct MultiPointZ{T,M} <: GeoInterface.AbstractMultiPoint
    MBR::Rect{T}
    points::Vector{Point{T}}
    zvalues::Vector{T}
    measures::Vector{M}
end

struct MultiPatch{T,M} <: GeoInterface.AbstractGeometry
    MBR::Rect{T}
    parts::Vector{Int32}
    parttypes::Vector{Int32}
    points::Vector{Point{T}}
    zvalues::Vector{T}
    # measures::Vector{M} # (optional)
end

const PointM{T, M} = MetaPoint{2, T, (:measure,), Tuple{M}}
const PointZ{T, M} = MetaPoint{3, T, (:measure,), Tuple{M}}


const SHAPETYPE = Dict{Int32, Any}(
    0  => NullShape,
    1  => Point{Float64},
    3  => Polyline{Float64},
    5  => Polygon{Float64},
    8  => MultiPoint{Float64},
    11 => PointZ{Float64,Float64},
    13 => PolylineZ{Float64,Float64},
    15 => PolygonZ{Float64,Float64},
    18 => MultiPointZ{Float64,Float64},
    21 => PointM{Float64, Float64},
    23 => PolylineM{Float64,Float64},
    25 => PolygonM{Float64,Float64},
    28 => MultiPointM{Float64,Float64},
    31 => MultiPatch{Float64,Float64}
)

struct Handle{T <: GeoInterface.AbstractGeometry}
    code::Int32
    length::Int32
    version::Int32
    shapeType::Int32
    MBR::Rect2D{Float64}
    zrange::Interval{Float64}
    mrange::Interval{Float64}
    shapes::Vector{T}
end
function read_t(io, ::Type{T}) where T
    ref = Ref{T}()
    read!(io, ref)
    ref[]
end
function Base.read(io::IO, T::Type{<: ShapeRect})
    rect = read_t(io, T)
    xy = (rect.left, rect.bottom)
    xymax = (rect.right, rect.top)
    Rect2D(xy..., (xymax .- xy)...)
end

Base.read(io::IO,::Type{NullShape}) = NullShape()


function read(io, ::Type{LineString{2, T}})
    numparts = read(io, Int32)
    numpoints = read(io, Int32)
    parts = Vector{LineFace{Int32}}(undef, numparts)
    read!(io, parts)
    points = Vector{Point{2, T}}(undef, numpoints)
    read!(io, points)
    return LineString(view(points, parts))
end

function Base.read(io::IO, PL::Type{LineString{T}}) where T
    box = read(io, Rect{T})
    LineString{T}(read_line!(io, PL))
end

function Base.read(io::IO,::Type{PolylineM{T,M}}) where {T,M}
    box = read(io, Rect{T})
    numparts = read(io,Int32)
    numpoints = read(io,Int32)
    parts = Vector{LineFace{Int32}}(undef, numparts)
    read!(io, parts)
    points = Vector{Point{2, T}}(undef, numpoints)
    read!(io, points)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    PolylineM{T,M}(faceview(structview(PointM{T, M}, points, measures), parts))
end

function Base.read(io::IO,::Type{PolylineZ{T,M}}) where {T,M}
    box = read(io,Rect{T})
    numparts = read(io,Int32)
    numpoints = read(io,Int32)
    parts = Vector{Int32}(undef, numparts)
    read!(io, parts)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    zrange = Vector{T}(undef, 2)
    read!(io, zrange)
    zvalues = Vector{T}(undef, numpoints)
    read!(io, zvalues)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    PolylineZ{T,M}(box,parts,points,zvalues,measures)
end

function Base.read(io::IO,::Type{Polygon{T}}) where T
    box = read(io,Rect{Float64})
    line = read_line!(io, LineString{2, T})
    Polygon{T}(box, line)
end

function Base.read(io::IO,::Type{PolygonM{T,M}}) where {T,M}
    box = read(io,Rect{Float64})
    numparts = read(io,Int32)
    numpoints = read(io,Int32)
    parts = Vector{Int32}(undef, numparts)
    read!(io, parts)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    PolygonM{T,M}(box,parts,points,measures)
end

function Base.read(io::IO,::Type{PolygonZ{T,M}}) where {T,M}
    box = read(io,Rect{Float64})
    numparts = read(io,Int32)
    numpoints = read(io,Int32)
    parts = Vector{Int32}(undef, numparts)
    read!(io, parts)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    zrange = Vector{T}(undef, 2)
    read!(io, zrange)
    zvalues = Vector{T}(undef, numpoints)
    read!(io, zvalues)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    PolygonZ{T,M}(box,parts,points,zvalues,measures)
end

function Base.read(io::IO,::Type{MultiPoint{T}}) where T
    box = read(io,Rect{Float64})
    numpoints = read(io,Int32)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    MultiPoint{T}(box,points)
end

function Base.read(io::IO,::Type{MultiPointM{T,M}}) where {T,M}
    box = read(io,Rect{Float64})
    numpoints = read(io,Int32)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    MultiPointM{T,M}(box,points,measures)
end

function Base.read(io::IO,::Type{MultiPointZ{T,M}}) where {T,M}
    box = read(io,Rect{Float64})
    numpoints = read(io,Int32)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    zrange = Vector{T}(undef, 2)
    read!(io, zrange)
    zvalues = Vector{T}(undef, numpoints)
    read!(io, zvalues)
    mrange = Vector{M}(undef, 2)
    read!(io, mrange)
    measures = Vector{M}(undef, numpoints)
    read!(io, measures)
    MultiPointZ{T,M}(box,points,zvalues,measures)
end

function Base.read(io::IO,::Type{MultiPatch{T,M}}) where {T,M}
    box = read(io,Rect{Float64})
    numparts = read(io,Int32)
    numpoints = read(io,Int32)
    parts = Vector{Int32}(undef, numparts)
    read!(io, parts)
    parttypes = Vector{Int32}(undef, numparts)
    read!(io, parttypes)
    points = Vector{Point{T}}(undef, numpoints)
    read!(io, points)
    zrange = Vector{T}(undef, 2)
    read!(io, zrange)
    zvalues = Vector{T}(undef, numpoints)
    read!(io, zvalues)
    # mrange = Vector{M}(2)
    # read!(io, mrange)
    # measures = Vector{M}(numpoints)
    # read!(io, measures)
    MultiPatch{T,M}(box,parts,parttypes,points,zvalues) #,measures)
end

function Base.read(io::IO,::Type{Handle})
    code = bswap(read(io,Int32))
    read!(io, Vector{Int32}(undef, 5))
    fileSize = bswap(read(io,Int32))
    version = read(io,Int32)
    shapeType = read(io,Int32)
    MBR = read(io,Rect{Float64})
    zmin = read(io,Float64)
    zmax = read(io,Float64)
    mmin = read(io,Float64)
    mmax = read(io,Float64)
    jltype = SHAPETYPE[shapeType]
    shapes = Vector{jltype}(undef, 0)
    file = Handle(code,fileSize,version,shapeType,MBR,Interval(zmin,zmax),Interval(mmin,mmax),shapes)
    while(!eof(io))
        num = bswap(read(io,Int32))
        rlength = bswap(read(io,Int32))
        shapeType = read(io,Int32)
        push!(shapes, read(io, jltype))
    end
    file
end
