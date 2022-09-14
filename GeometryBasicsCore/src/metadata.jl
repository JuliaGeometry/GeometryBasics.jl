"""
    MetaType(::Type{T})

Returns the Meta Type corresponding to `T`
E.g:
```julia
MetaType(Point) == PointMeta
```
"""
MetaType(::Type{T}) where {T} = error("No Meta Type for $T")

"""
    MetaFree(::Type{T})

Returns the original type containing no metadata for `T`
E.g:
```julia
MetaFree(PointMeta) == Point
```
"""
MetaFree(::Type{T}) where {T} = error("No meta free Type for $T")

"""
    meta(x::MetaObject)

Returns the metadata of `x`
"""
meta(x::T) where {T} = error("$T has no meta!")

metafree(x::T) where {T} = x


"""

    MetaT(geometry, meta::NamedTuple)
    MetaT(geometry; meta...)
    
Returns a `MetaT` that holds a geometry and its metadata

`MetaT` acts the same as `Meta` method.
The difference lies in the fact that it is designed to handle
geometries and metadata of different/heterogeneous types.

eg: While a Point MetaGeometry is a `PointMeta`, the MetaT representation is `MetaT{Point}`
The downside being it's not subtyped to `AbstractPoint` like a `PointMeta` is.

Example:
```julia
julia> MetaT(Point(1, 2), city = "Mumbai")
MetaT{Point{2,Int64},(:city,),Tuple{String}}([1, 2], (city = "Mumbai",))
```
"""
struct MetaT{T,Names,Types}
    main::T
    meta::NamedTuple{Names,Types}
end

MetaT(x; kwargs...) = MetaT(x, values(kwargs))

"""

    metafree(x::MetaT)
    metafree(x::Array{MetaT})

Free the MetaT from metadata
i.e. returns the geometry/array of geometries
"""
function metafree(x::MetaT)
    return getfield(x, :main)
end
metafree(x::AbstractVector{<:MetaT}) = map(metafree, x)

"""

    meta(x::MetaT)
    meta(x::Array{MetaT})
    
Returns the metadata of a `MetaT`
"""
function meta(x::MetaT)
    return getfield(x, :meta)
end
meta(x::AbstractVector{<:MetaT}) = map(meta, x)

# helper methods
function Base.getproperty(x::MetaT, field::Symbol)
    return if field == :main
        metafree(x)
    elseif field == :meta
        meta(x)
    else
        getproperty(meta(x), field)
    end
end

Base.propertynames(x::MetaT) = (:main, propertynames(meta(x))...)
getnamestypes(::Type{MetaT{T,Names,Types}}) where {T,Names,Types} = (T, Names, Types)

Base.getindex(x::MetaT, idx::Int) = getindex(metafree(x), idx)
Base.size(x::MetaT) = size(metafree(x))
