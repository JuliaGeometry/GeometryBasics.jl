##
# Generic base overloads
Base.extrema(primitive::GeometryPrimitive) = (minimum(primitive), maximum(primitive))
function widths(x::AbstractRange)
    mini, maxi = Float32.(extrema(x))
    return maxi - mini
end

##
# conversion & decompose
convert_simplex(::Type{T}, x::T) where {T} = (x,)
convert_simplex(::Type{Vec{N, T}}, x::Vec{N, T}) where {N, T} = x

function convert_simplex(NFT::Type{NgonFace{N,T1}},
                         f::Union{NgonFace{N,T2}}) where {T1,T2,N}
    return (convert(NFT, f),)
end

function convert_simplex(NFT::Type{NgonFace{3,T}}, f::NgonFace{3,T2}) where {T,T2}
    return (convert(NFT, f),)
end
function convert_simplex(NFT::Type{NgonFace{2,T}}, f::NgonFace{2,T2}) where {T,T2}
    return (convert(NFT, f),)
end

"""
    convert_simplex(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function convert_simplex(::Type{TriangleFace{T}},
                                    f::Union{SimplexFace{N},NgonFace{N}}) where {T,N}
    3 <= N || error("decompose not implemented for N <= 3 yet. N: $N")# other wise degenerate
    v = Expr(:tuple)
    for i in 3:N
        push!(v.args, :(TriangleFace{T}(f[1], f[$(i - 1)], f[$i])))
    end
    return v
end

"""
    convert_simplex(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function convert_simplex(::Type{LineFace{T}},
                                    f::Union{SimplexFace{N},NgonFace{N}}) where {T,N}
    2 <= N || error("decompose not implented for N <= 2 yet. N: $N")# other wise degenerate

    v = Expr(:tuple)
    for i in 1:(N - 1)
        push!(v.args, :(LineFace{$T}(f[$i], f[$(i + 1)])))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace{$T}(f[$N], f[1])))
    return v
end

to_pointn(::Type{T}, x) where {T<:Point} = convert_simplex(T, x)[1]

# disambiguation method overlords
convert_simplex(::Type{Point}, x::Point) = (x,)
convert_simplex(::Type{Point{N,T}}, p::Point{N,T}) where {N,T} = (p,)
function convert_simplex(::Type{Point{N,T}}, x) where {N,T}
    N2 = length(x)
    return (Point{N,T}(ntuple(i -> i <= N2 ? T(x[i]) : T(0), N)),)
end

function convert_simplex(::Type{Vec{N,T}}, x) where {N,T}
    N2 = length(x)
    return (Vec{N,T}(ntuple(i -> i <= N2 ? T(x[i]) : T(0), N)),)
end

collect_with_eltype(::Type{T}, vec::Vector{T}) where {T} = vec
collect_with_eltype(::Type{T}, vec::AbstractVector{T}) where {T} = collect(vec)

function collect_with_eltype(::Type{T}, iter) where {T}
    isempty(iter) && return T[]
    # We need to get `eltype` information from `iter`, it seems to be `Any`
    # most of the time so the eltype checks here don't actually work
    l = if Base.IteratorSize(iter) isa Union{Base.HasShape,Base.HasLength}
        if Base.IteratorEltype(iter) isa Base.HasEltype && isconcretetype(eltype(iter))
            # Work out the exact length
            length(convert_simplex(T, first(iter))) * length(iter)
        else
            # We know it is at least the length of iter,
            # after that we will `push!` if we have to
            length(iter)
        end
    else
        0
    end
    n = 0
    result = Vector{T}(undef, l)
    for element in iter
        # convert_simplex always returns a tuple,
        # so that e.g. convert(Triangle, quad) can return 2 elements
        for telement in convert_simplex(T, element)
            n += 1
            if n > l
                push!(result, telement)
            else
                result[n] = telement
            end
        end
    end
    return result
end

"""
The unnormalized normal of three vertices.
"""
function orthogonal_vector(v1, v2, v3)
    a = v2 .- v1
    b = v3 .- v1
    return cross(a, b)
end

"""
```
normals{VT,FD,FT,FO}(vertices::Vector{Point{3, VT}},
                    faces::Vector{Face{FD,FT,FO}},
                    NT = Normal{3, VT})
```
Compute all vertex normals.
"""
function normals(vertices::AbstractVector{Point{3,T}}, faces::AbstractVector{F};
                 normaltype=Vec{3,T}) where {T,F<:NgonFace}
    return normals(vertices, faces, normaltype)
end

function normals(vertices::AbstractVector{<:Point{3}}, faces::AbstractVector{F},
                 ::Type{NormalType}) where {F<:NgonFace,NormalType}
    normals_result = zeros(NormalType, length(vertices))
    for face in faces
        v = vertices[face]
        # we can get away with two edges since faces are planar.
        n = orthogonal_vector(v[1], v[2], v[3])
        for i in 1:length(F)
            fi = face[i]
            normals_result[fi] = normals_result[fi] .+ n
        end
    end
    normals_result .= normalize.(normals_result)
    return normals_result
end
