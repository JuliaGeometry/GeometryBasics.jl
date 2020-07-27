struct metag{T, Names, Types}
    data::T
    rest::NamedTuple{Names, Types}
end

metag(x; kwargs...) = metag(x, values(kwargs))

Base.getproperty(mg::metag, s::Symbol) = s == :data ? getfield(mg, 1) : getproperty(getfield(mg, 2), s)

getnamestypes(::Type{metag{T, Names, Types}}) where {T, Names, Types} = (T, Names, Types)

function StructArrays.staticschema(::Type{mg}) where {mg<:metag}
    T, names, types = getnamestypes(mg)
    NamedTuple{(:data, names...), Base.tuple_type_cons(T, types)}
end

function StructArrays.createinstance(::Type{mg}, x, args...) where {mg<:metag}
    T, names, types = getnamestypes(mg)
    metag(x, NamedTuple{names, types}(args))
end
