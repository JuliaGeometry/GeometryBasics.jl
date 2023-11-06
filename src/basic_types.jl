"""
Abstract Geometry in R{Dim} with Number type T
"""
abstract type AbstractGeometry{Dim,T<:Number} end
Base.ndims(::AbstractGeometry{Dim}) where {Dim} = Dim

"""
Atomic geometry primitive in R{Dim}  with Number type T
"""
abstract type GeometryPrimitive{Dim,T} <: AbstractGeometry{Dim,T} end

"""
    AbstractMesh

An abstract mesh is a collection of `GeometryPrimitive` elements (e.g. simplices / ngons).
The connections are defined via faces(mesh), the coordinates of the elements are returned by
coordinates(mesh). Arbitrary meta information can be attached per point or per face.
"""
abstract type AbstractMesh{Dim, T} <: AbstractGeometry{Dim, T} end

abstract type AbstractFace{N,T} <: StaticVector{N,T} end
