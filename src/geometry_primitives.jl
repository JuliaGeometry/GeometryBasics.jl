##
# Generic base overloads
Base.extrema(primitive::GeometryPrimitive) = (minimum(primitive), maximum(primitive))
function widths(x::AbstractRange)
    mini, maxi = Float32.(extrema(x))
    return maxi - mini
end

##
# conversion & decompose
convert_simplex(::Type{T}, x::T) where T = (x,)

function convert_simplex(NFT::Type{NgonFace{N, T1}}, f::Union{NgonFace{N, T2}}) where {T1, T2, N}
    return (convert(NFT, f),)
end

convert_simplex(NFT::Type{NgonFace{3,T}}, f::NgonFace{3,T2}) where {T, T2} = (convert(NFT, f),)
convert_simplex(NFT::Type{NgonFace{2,T}}, f::NgonFace{2,T2}) where {T, T2} = (convert(NFT, f),)

"""
    convert_simplex(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function convert_simplex(::Type{TriangleFace{T}}, f::Union{SimplexFace{N}, NgonFace{N}}) where {T, N}
    3 <= N || error("decompose not implemented for N <= 3 yet. N: $N")# other wise degenerate
    v = Expr(:tuple)
    for i = 3:N
        push!(v.args, :(TriangleFace{T}(f[1], f[$(i-1)], f[$i])))
    end
    return v
end

"""
    convert_simplex(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function convert_simplex(::Type{LineFace{T}}, f::Union{SimplexFace{N}, NgonFace{N}}) where {T, N}
    2 <= N || error("decompose not implented for N <= 2 yet. N: $N")# other wise degenerate

    v = Expr(:tuple)
    for i = 1:N-1
        push!(v.args, :(LineFace{$T}(f[$i], f[$(i+1)])))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace{$T}(f[$N], f[1])))
    return v
end

to_pointn(::Type{T}, x) where T<:Point = convert_simplex(T, x)[1]

# disambiguation method overlords
convert_simplex(::Type{Point}, x::Point) = (x,)
convert_simplex(::Type{Point{N,T}}, p::Point{N,T}) where {N, T} = (p,)
function convert_simplex(::Type{Point{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Point{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end

function convert_simplex(::Type{Vec{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Vec{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end

collect_with_eltype(::Type{T}, vec::Vector{T}) where T = vec
collect_with_eltype(::Type{T}, vec::AbstractVector{T}) where T = collect(vec)

function collect_with_eltype(::Type{T}, iter) where T
    # TODO we could be super smart about allocating the right length
    # but its kinda annoying, since e.g. T == Triangle and first(iter) isa Quad
    # will need double the length etc - but could all be figured out ;)
    result = T[]
    for element in iter
        # convert_simplex always returns a tuple,
        # so that e.g. convert(Triangle, quad) can return 2 elements
        for telement in convert_simplex(T, element)
            push!(result, telement)
        end
    end
    return result
end


"""
The unnormalized normal of three vertices.
"""
function orthogonal_vector(v1, v2, v3)
    a = v2 - v1
    b = v3 - v1
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
function normals(vertices::AbstractVector{<: AbstractPoint{3, T}},
                 faces::AbstractVector{F};
                 normaltype=Vec{3, T}) where {T, F <: NgonFace}
    normals_result = zeros(normaltype, length(vertices)) # initilize with same type as verts but with 0
    for face in faces
        v = metafree.(vertices[face])
        # we can get away with two edges since faces are planar.
        n = orthogonal_vector(v[1], v[2], v[3])
        for i =1:length(F)
            fi = face[i]
            normals_result[fi] = normals_result[fi] + n
        end
    end
    normals_result .= normalize.(normals_result)
    return normals_result
end
