##
# Generic base overloads

Base.extrema(primitive::GeometryPrimitive) = (minimum(primitive), maximum(primitive))
function widths(x::AbstractRange)
    mini, maxi = Float32.(extrema(x))
    maxi - mini
end

##
# conversion & decompose

"""
    simplex_convert(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function simplex_convert(::Type{TriangleFace{T}}, f::NgonFace{N}) where {T, N}
    3 <= N || error("decompose not implented for N <= 3 yet. N: $N")# other wise degenerate
    v = Expr(:tuple)
    for i = 3:N
        push!(v.args, :(TriangleFace{T}(f[1], f[$(i-1)], f[$i])))
    end
    return v
end

"""
    simplex_convert(::Type{Face{2}}, f::Face{N})

Extract all line segments in a Face.
"""
@generated function simplex_convert(::Type{LineFace{T}}, f::NgonFace{N}) where {T, N}
    2 <= N || error("decompose not implented for N <= 2 yet. N: $N")# other wise degenerate

    v = Expr(:tuple)
    for i = 1:N-1
        push!(v.args, :(LineFace{$T}(f[$i], f[$(i+1)])))
    end
    # connect vertices N and 1
    push!(v.args, :(LineFace{$T}(f[$N], f[1])))
    return v
end

function simplex_convert(::Type{Point{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Point{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end

function simplex_convert(::Type{Vec{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Vec{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end


function collect_with_eltype(::Type{T}, iter) where T
    # TODO we could be super smart about allocating the right length
    # but its kinda annoying, since e.g. T == Triangle and first(iter) isa Quad
    # will need double the length etc - but could all be figured out ;)
    result = T[]
    for element in iter
        # simplex_convert always returns a tuple,
        # so that e.g. convert(Triangle, quad) can return 2 elements
        for telement in simplex_convert(T, element)
            push!(result, telement)
        end
    end
    return result
end

function decompose(::Type{T}, primitive::AbstractVector{T}) where {T}
    return primitive
end

function decompose(::Type{P}, primitive) where {P<:AbstractPoint}
    return collect_with_eltype(P, coordinates(primitive))
end

function decompose(::Type{F}, primitive) where {F<:AbstractFace}
    return collect_with_eltype(F, faces(primitive))
end

function decompose_uv(primitive::GeometryPrimitive)
    return collect_with_eltype(Vec2f0, texturecoordinates(primitive))
end

function decompose_normals(primitive::GeometryPrimitive)
    return collect_with_eltype(Vec3f0, normals(primitive))
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
function normals(vertices::AbstractVector{Point{3, T}}, faces::AbstractVector{F}) where {T, F <: NgonFace}
    normals_result = zeros(Vec{3, T}, length(vertices)) # initilize with same type as verts but with 0
    for face in faces
        v = vertices[face]
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


##
# Some more primitive types


"""
A `HyperSphere` is a generalization of a sphere into N-dimensions.
A `center` and radius, `r`, must be specified.
"""
struct HyperSphere{N, T} <: GeometryPrimitive{N, T}
    center::Point{N, T}
    r::T
end
"""
An alias for a HyperSphere of dimension 2. i.e. `Circle{T}` -> `HyperSphere{2, T}`
"""
const Circle{T} = HyperSphere{2, T}

"""
An alias for a HyperSphere of dimension 3. i.e. `Sphere{T}` -> `HyperSphere{3, T}`
"""
const Sphere{T} = HyperSphere{3, T}

"""
A rectangle in 3D space.
"""
struct Quad{T} <: GeometryPrimitive{3, T}
    downleft::Vec{3, T}
    width   ::Vec{3, T}
    height  ::Vec{3, T}
end

struct Pyramid{T} <: GeometryPrimitive{3, T}
    middle::Point{3, T}
    length::T
    width ::T
end

struct Particle{N, T} <: GeometryPrimitive{N, T}
    position::Point{N, T}
    velocity::Vec{N, T}
end

"""
A `Cylinder` is a 2D rectangle or a 3D cylinder defined by its origin point,
its extremity and a radius. `origin`, `extremity` and `r`, must be specified.
"""
struct Cylinder{N, T} <: GeometryPrimitive{N, T}
    origin::Point{N,T}
    extremity::Point{N,T}
    r::T
end

"""
A `Cylinder2` or `Cylinder3` is a 2D/3D cylinder defined by its origin point,
its extremity and a radius. `origin`, `extremity` and `r`, must be specified.
"""
const Cylinder2{T} = Cylinder{2, T}
const Cylinder3{T} = Cylinder{3, T}


origin(c::Cylinder{N, T}) where {N, T} = c.origin
extremity(c::Cylinder{N, T}) where {N, T} = c.extremity
radius(c::Cylinder{N, T}) where {N, T} = c.r
height(c::Cylinder{N, T}) where {N, T} = norm(c.extremity - c.origin)
direction(c::Cylinder{N, T}) where {N, T} = (c.extremity .- c.origin) ./ height(c)

function rotation(c::Cylinder{2, T}) where T
    d2 = direction(c); u = @SVector [d2[1], d2[2], T(0)]
    v = @MVector [u[2], -u[1], T(0)]
    normalize!(v)
    return hcat(v, u, @SVector T[0, 0, 1])
end

function rotation(c::Cylinder{3, T}) where T
    d3 = direction(c); u = @SVector [d3[1], d3[2], d3[3]]
    if abs(u[1]) > 0 || abs(u[2]) > 0
        v = @MVector [u[2], -u[1], T(0)]
    else
        v = @MVector [T(0), -u[3], u[2]]
    end
    normalize!(v)
    w = @SVector [u[2] * v[3] - u[3] * v[2], -u[1] * v[3] + u[3] * v[1], u[1] * v[2] - u[2] * v[1]]
    return hcat(v, w, u)
end

function coordinates(s::Cylinder{2}, nvertices=24)
    r = Rect(c.origin[1] - c.r/2, c.origin[2], c.r, height(c))
    M = rotation(c)
    points = coordinates(r, nvertices)
    vo = origin(s)
    return (M * ((point .- vo) .+ vo) for point in points)
end

function coordinates(c::Cylinder{3, T}, nvertices=30) where T
    isodd(nvertices) && (nvertices = 2 * div(nvertices, 2))
    nbv = div(nvertices, 2)
    M = rotation(c)
    h = height(c)
    range = 1:(2 * nbv + 2)
    function inner(i)
        phi = T((2π * (j - 1)) / nbv)
        if i == length(range)
            return c.extremity
        elseif i == length(range) - 1
            origin(c)
        elseif isodd(i)
            return PT(M * PT(c.r * cos(phi), c.r * sin(phi),0)) .+ PT(c.origin)
        else
            PT(M * PT(c.r * cos(phi), c.r * sin(phi),h)) .+ PT(c.origin)
        end
    end

    return (inner(i) for i in range)
end


function coordinates(c::Cylinder{3}, nvertices=30)
    isodd(nvertices) ? nvertices = 2 * div(nvertices, 2) : nothing
    nvertices < 8 ? nvertices = 8 : nothing; nbv = Int(nvertices / 2)
    indexes = Vector{GLTriangle}(undef, nvertices)
    index = 1
    for j = 1:(nbv-1)
        indexes[index] = (index + 2, index + 1, index)
        indexes[index + 1] = ( index + 3, index + 1, index + 2)
        index += 2
    end
    indexes[index] = (1, index + 1, index)
    indexes[index + 1] = (2, index + 1, 1)
    for i in 1:length(indexes)
        if i%2 == 1
            push!(indexes, (indexes[i][1], indexes[i][3], 2*nbv+1))
        else
            push!(indexes, (indexes[i][2], indexes[i][1], 2*nbv+2))
        end
    end
    return indexes
end


function texturecoordinates(s::Cylinder, nvertices=24)
    ux = LinRange(0, 1, nvertices)
    return ivec((Vec(φ, θ) for θ in reverse(ux), φ in ux))
end

function faces(sphere::Cylinder, nvertices=24)
    return faces(Rect(0, 0, 1, 1), (nvertices, nvertices))
end

function normals(s::Cylinder{T}, nvertices=24) where T
    return coordinates(Sphere(Point{3, T}(0), 1f0), nvertices)
end


##
# Sphere

HyperSphere{N}(p::Point{N, T}, number) where {N, T} = HyperSphere{N, T}(p, convert(T, number))

widths(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(radius(c)*2)
radius(c::HyperSphere) = c.r
origin(c::HyperSphere) = c.center
minimum(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(origin(c)) - Vec{N, T}(radius(c))
maximum(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(origin(c)) + Vec{N, T}(radius(c))
function isinside(c::Circle, x::Real, y::Real)
    @inbounds ox, oy = origin(c)
    xD = abs(ox - x)
    yD = abs(oy - y)
    return xD <= c.r && yD <= c.r
end

centered(S::Type{HyperSphere{N, T}}) where {N, T} = S(Vec{N,T}(0), T(0.5))
centered(::Type{T}) where {T <: HyperSphere} = centered(HyperSphere{ndims_or(T, 3), eltype_or(T, Float32)})

function coordinates(s::Circle, nvertices=64)
    rad = radius(s)
    inner(fi) = Point(rad*sin(fi + pi), rad*cos(fi + pi)) .+ origin(s)
    return (inner(fi) for fi in LinRange(0, 2pi, nvertices))
end

function coordinates(s::Sphere, n = 24, nvertices=24)
    θ = LinRange(0, pi, nvertices); φ = LinRange(0, 2pi, nvertices)
    inner(θ, φ) = Point(cos(φ)*sin(θ), sin(φ)*sin(θ), cos(θ)) .* s.r .+ s.center
    return ivec((inner(θ, φ) for θ in θ, φ in φ))
end


function texturecoordinates(s::Sphere, nvertices=24)
    ux = LinRange(0, 1, nvertices)
    return ivec(((φ, θ) for θ in reverse(ux), φ in ux))
end

function faces(sphere::Sphere, nvertices=24)
    return faces(Rect(0, 0, 1, 1), (nvertices, nvertices))
end

function normals(s::Sphere{T}, nvertices=24) where {T}
    return coordinates(Sphere(Point{3, T}(0), 1), nvertices)
end
