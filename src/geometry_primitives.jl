##
# Generic base overloads
Base.extrema(primitive::GeometryPrimitive) = (minimum(primitive), maximum(primitive))
function widths(x::AbstractRange)
    mini, maxi = Float32.(extrema(x))
    return maxi - mini
end

##
# conversion & decompose
"""
    convert_simplex(::Type{TargetType}, x)

Used to convert one object into another in `decompose(::Type{TargetType}, xs)`.
"""
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

# TODO: generic?
function convert_simplex(::Type{TriangleFace{T}}, f::SimplexFace{4}) where {T}
    TF = TriangleFace{T}
    return (TF(f[2],f[3],f[4]), TF(f[1],f[3],f[4]), TF(f[1],f[2],f[4]), TF(f[1],f[2],f[3]))
end

"""
    convert_simplex(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function convert_simplex(::Type{LineFace{T}}, f::NgonFace{N}) where {T,N}
    2 <= N || error("decompose not implemented for N <= 2 yet. N: $N")# other wise degenerate

    v = Expr(:tuple)
    for i in 1:(N - 1)
        push!(v.args, :(LineFace{$T}(f[$i], f[$(i + 1)])))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace{$T}(f[$N], f[1])))
    return v
end

@generated function convert_simplex(::Type{LineFace{T}}, f::SimplexFace{N}) where {T,N}
    2 <= N || error("decompose not implemented for N <= 2 yet. N: $N")# other wise degenerate

    v = Expr(:tuple)
    for i in 1:(N - 1)
        for j in i+1:N
            push!(v.args, :(LineFace{$T}(f[$i], f[$j])))
        end
    end
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
collect_with_eltype(::Type{T}, vec::FaceView{T}) where {T} = vec

function collect_with_eltype(::Type{T}, iter) where {T}
    return collect_with_eltype!(Vector{T}(undef, 0), iter)
end
function collect_with_eltype(::Type{T}, iter::FaceView) where {T}
    return FaceView(collect_with_eltype!(Vector{T}(undef, 0), iter.data), iter.faces)
end

function collect_with_eltype!(target::AbstractVector{T}, vec::AbstractVector{T}) where {T}
    return append!(target, vec)
end

function collect_with_eltype!(result::AbstractVector{T}, iter) where {T}
    isempty(iter) && return result

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

    # Allow result to be pre-filled for handling faces with mesh.views
    sizehint!(result, length(result) + l)

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
    orthogonal_vector(vertices)

Calculates an orthogonal vector to a face or polygon defined by a set of 
ordered points contained in the point vector `vertices`. Note that for 2D points
a scalar is returned (the dot product of the orthogonal vector with the 
positive Z-direction). 
"""
function orthogonal_vector(::Type{VT},vertices) where VT     
    c = zeros(VT) # Inherit vector type from input 
    j = length(vertices)
    @inbounds for i in eachindex(vertices) # Use shoelace approach
        c += cross(vertices[j],vertices[i]) # Add each edge contribution          
        j = i
    end
    return c
end

function orthogonal_vector(vertices)
    return orthogonal_vector(eltype(vertices),vertices)
end

"""
    normals(positions::Vector{Point3{T}}, faces::Vector{<: NgonFace}[; normaltype = Vec3{T}])

Compute vertex normals from the given `positions` and `faces`.

This runs through all faces, computing a face normal each and adding it to every
involved vertex. The direction of the face normal is based on winding direction 
and assumed counter-clockwise faces. At the end the summed face normals are 
normalized again to produce a vertex normal.
"""
function normals(vertices::AbstractVector{Point{3,T}}, faces::AbstractVector{F};
                 normaltype=Vec{3,T}) where {T,F<:NgonFace}
    return normals(vertices, faces, normaltype)
end

function normals(vertices::AbstractVector{<:Point{3}}, faces::AbstractVector{<: NgonFace},
                 ::Type{NormalType}) where {NormalType}
                                  
    normals_result = zeros(NormalType, length(vertices))
    for face in faces
        v = coordinates(vertices[face])
        # we can get away with two edges since faces are planar.
        n = orthogonal_vector(NormalType, v)
        for i in 1:length(face)
            fi = face[i]
            normals_result[fi] = normals_result[fi] .+ n
        end
    end
    normals_result .= normalize.(normals_result)
    return normals_result
end


"""
    face_normals(positions::Vector{Point3{T}}, faces::Vector{<: NgonFace}[, target_type = Vec3{T}])

Compute face normals from the given `positions` and `faces` and returns an
appropriate `FaceView`.
"""
function face_normals(
        positions::AbstractVector{<:Point3{T}}, fs::AbstractVector{<: AbstractFace}; 
        normaltype = Vec3{T}) where {T}
    return face_normals(positions, fs, normaltype)
end
    
@generated function face_normals(positions::AbstractVector{<:Point3}, fs::AbstractVector{F},
        ::Type{NormalType}) where {F<:NgonFace,NormalType}
    
    # If the facetype is not concrete it likely varies and we need to query it 
    # doing the iteration
    FT = ifelse(isconcretetype(F), :($F), :(typeof(f)))

    quote
        normals = resize!(NormalType[], length(fs))
        faces   = resize!(F[], length(fs))

        for (i, f) in enumerate(fs)
            ps = coordinates(positions[f])
            n = orthogonal_vector(NormalType, ps)
            normals[i] = normalize(n)
            faces[i] = $(FT)(i)
        end

        return FaceView(normals, faces)
    end
end
