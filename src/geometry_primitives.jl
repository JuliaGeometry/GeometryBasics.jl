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
    convert_simplex(::Type{Face{3}}, f::Face{N})

Triangulate an N-Face into a tuple of triangular faces.
"""
@generated function convert_simplex(::Type{TriangleFace{T}}, f::NgonFace{N}) where {T, N}
    3 <= N || error("decompose not implented for N <= 3 yet. N: $N")# other wise degenerate
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
@generated function convert_simplex(::Type{LineFace{T}}, f::NgonFace{N}) where {T, N}
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

convert_simplex(::Type{Point}, x::Point) = (x,)
function convert_simplex(::Type{Point{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Point{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end

function convert_simplex(::Type{Vec{N, T}}, x) where {N, T}
    N2 = length(x)
    return (Vec{N, T}(ntuple(i-> i <= N2 ? T(x[i]) : T(0), N)),)
end

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


function faces(primitive, nvertices=30)
    # doesn't have any specific algorithm to generate faces
    # so will try to triangulate the coordinates!
    return nothing
end

function decompose(::Type{F}, primitive, args...) where {F<:AbstractFace}
    f = faces(primitive, args...)
    f === nothing && return nothing
    return collect_with_eltype(F, f)
end

function decompose(::Type{T}, primitive::AbstractVector{T}) where {T}
    return primitive
end

function decompose(::Type{T}, primitive::AbstractVector{T}) where {T<:AbstractFace}
    return primitive
end

function decompose(::Type{T}, primitive::AbstractVector{T}) where {T<:AbstractPoint}
    return primitive
end

function decompose(::Type{P}, primitive, args...) where {P<:AbstractPoint}
    return collect_with_eltype(P, coordinates(primitive, args...))
end

function decompose(::Type{Point}, primitive::Union{GeometryPrimitive{Dim}, Mesh{Dim}}, args...) where {Dim}
    return collect_with_eltype(Point{Dim, Float32}, coordinates(primitive, args...))
end

# Dispatch type to make `decompose(UV{Vec2f0}, priomitive)` work
struct UV{T} end
UV(::Type{T}) where T = UV{T}()
struct UVW{T} end
UVW(::Type{T}) where T = UVW{T}()
struct Normal{T} end
Normal(::Type{T}) where T = Normal{T}()

function decompose(::UV{T}, primitive::GeometryPrimitive, args...) where T
    return collect_with_eltype(T, texturecoordinates(primitive, args...))
end

decompose_uv(args...) = decompose(UV(Vec2f0), args...)

function decompose(::UVW{T}, primitive::GeometryPrimitive, args...) where T
    return collect_with_eltype(T, texturecoordinates(primitive, args...))
end

function normals(primitive, nvertices=30)
    # doesn't have any specific algorithm to generate normals
    # so will be generated from faces + positions
    return nothing
end


function decompose(::Normal{T}, primitive::GeometryPrimitive, args...) where T
    n = normals(primitive, args...)
    n === nothing && return nothing
    return collect_with_eltype(T, n)
end

decompose_normals(args...) = decompose(Normal(Vec3f0), args...)

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

function normals(mesh::AbstractMesh)
    if hasproperty(mesh, :normals)
        return mesh.normals
    else
        return normals(coordinates(mesh), faces(mesh))
    end
end

##
# Some more primitive types

"""
    HyperSphere{N, T}

A `HyperSphere` is a generalization of a sphere into N-dimensions.
A `center` and radius, `r`, must be specified.
"""
struct HyperSphere{N, T} <: GeometryPrimitive{N, T}
    center::Point{N, T}
    r::T
end
"""
    Circle{T}

An alias for a HyperSphere of dimension 2. (i.e. `HyperSphere{2, T}`)
"""
const Circle{T} = HyperSphere{2, T}

"""
    Sphere{T}

An alias for a HyperSphere of dimension 3. (i.e. `HyperSphere{3, T}`)
"""
const Sphere{T} = HyperSphere{3, T}

"""
    Quad{T}

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
    Cylinder{N, T}

A `Cylinder` is a 2D rectangle or a 3D cylinder defined by its origin point,
its extremity and a radius. `origin`, `extremity` and `r`, must be specified.
"""
struct Cylinder{N, T} <: GeometryPrimitive{N, T}
    origin::Point{N,T}
    extremity::Point{N,T}
    r::T
end

"""
    Cylinder2{T}
    Cylinder3{T}

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

best_nvertices(x::Cylinder{2}) = (2, 2)

function coordinates(c::Cylinder{2, T}, nvertices=(2, 2)) where T
    r = Rect(c.origin[1] - c.r/2, c.origin[2], c.r, height(c))
    M = rotation(c)
    points = coordinates(r, nvertices)
    vo = to_pointn(Point3{T}, origin(c))
    return (M * (to_pointn(Point3{T}, point) .- vo) .+ vo for point in points)
end

function faces(sphere::Cylinder{2}, nvertices=(2, 2))
    return faces(Rect(0, 0, 1, 1), nvertices)
end

best_nvertices(x::Cylinder{3}) = 30

function coordinates(c::Cylinder{3, T}, nvertices=30) where T
    if isodd(nvertices)
        nvertices = 2 * (nvertices ÷ 2)
    end
    nvertices = max(8, nvertices);
    nbv = nvertices ÷ 2

    M = rotation(c)
    h = height(c)
    range = 1:(2 * nbv + 2)
    function inner(i)
        if i == length(range)
            return c.extremity
        elseif i == length(range) - 1
            return origin(c)
        else
            phi = T((2π * (((i + 1) ÷ 2) - 1)) / nbv)
            up = ifelse(isodd(i), 0, h)
            return (M * Point(c.r * cos(phi), c.r * sin(phi), up)) .+ c.origin
        end
    end

    return (inner(i) for i in range)
end

function faces(c::Cylinder{3}, facets=30)
    isodd(facets) ? facets = 2 * div(facets, 2) : nothing
    facets < 8 ? facets = 8 : nothing; nbv = Int(facets / 2)
    indexes = Vector{TriangleFace{Int}}(undef, facets)
    index = 1
    for j = 1:(nbv-1)
        indexes[index] = (index + 2, index + 1, index)
        indexes[index + 1] = ( index + 3, index + 1, index + 2)
        index += 2
    end
    indexes[index] = (1, index + 1, index)
    indexes[index + 1] = (2, index + 1, 1)

    for i = 1:length(indexes)
        i%2 == 1 ? push!(indexes, (indexes[i][1], indexes[i][3], 2*nbv+1)) : push!(indexes,(indexes[i][2], indexes[i][1], 2*nbv+2))
    end
    return indexes
end

##
# Sphere

HyperSphere{N}(p::Point{N, T}, number) where {N, T} = HyperSphere{N, T}(p, convert(T, number))

widths(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(radius(c)*2)
radius(c::HyperSphere) = c.r
origin(c::HyperSphere) = c.center

Base.minimum(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(origin(c)) - Vec{N, T}(radius(c))
Base.maximum(c::HyperSphere{N, T}) where {N, T} = Vec{N, T}(origin(c)) + Vec{N, T}(radius(c))

function Base.in(x::AbstractPoint{2}, c::Circle)
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

function texturecoordinates(s::Circle, nvertices=64)
    return coordinates(Circle(Point2f0(0.5), 0.5f0), nvertices)
end

function coordinates(s::Sphere, nvertices=24)
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

function coordinates(p::Pyramid{T}, nvertices=nothing) where {T}
    leftup = Point{3, T}(-p.width , p.width, 0) / 2
    leftdown = Point(-p.width, -p.width, 0) / 2
    tip = Point{3, T}(p.middle + Point{3, T}(0, 0, p.length))
    lu = Point{3, T}(p.middle + leftup)
    ld = Point{3, T}(p.middle + leftdown)
    ru = Point{3, T}(p.middle - leftdown)
    rd = Point{3, T}(p.middle - leftup)
    return Point{3, T}[
        tip, rd, ru,
        tip, ru, lu,
        tip, lu, ld,
        tip, ld, rd,
        rd,  ru, lu,
        lu,  ld, rd
    ]
end

function faces(r::Pyramid, nvertices=nothing) where FT
    return (TriangleFace(triangle) for triangle in TupleView{3}(1:18))
end
