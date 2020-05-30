# Decomposition


## GeometryBasic Mesh interface

GeometryBasic defines an interface, to decompose abstract geometries into
points and triangle meshes.
This can be done for any arbitrary primitive, by overloading the following interface:

```julia

function GeometryBasics.coordinates(rect::Rect2D, nvertices=(2,2))
    mini, maxi = extrema(rect)
    xrange, yrange = LinRange.(mini, maxi, nvertices)
    return ivec(((x,y) for x in xrange, y in yrange))
end

function GeometryBasics.faces(rect::Rect2D, nvertices=(2, 2))
    w, h = nvertices
    idx = LinearIndices(nvertices)
    quad(i, j) = QuadFace{Int}(idx[i, j], idx[i+1, j], idx[i+1, j+1], idx[i, j+1])
    return ivec((quad(i, j) for i=1:(w-1), j=1:(h-1)))
end
```
Those methods, for performance reasons, expect you to return an iterator, to make
materializing them with different element types allocation free. But of course,
can also return any `AbstractArray`.

With these methods defined, this constructor will magically work:

```julia
rect = Rect2D(0.0, 0.0, 1.0, 1.0)
m = GeometryBasics.mesh(rect)
```
If you want to set the `nvertices` argument, you need to wrap your primitive in a `Tesselation`
object:
```julia
m = GeometryBasics.mesh(Tesselation(rect, (50, 50)))
length(coordinates(m)) == 50^2
```

As you can see, `coordinates` and `faces` is also defined on a mesh
```julia
coordinates(m)
faces(m)
```
But will actually not be an iterator anymore. Instead, the mesh constructor uses
the `decompose` function, that will collect the result of coordinates and will
convert it to a concrete element type:
```julia
decompose(Point2f0, rect) == convert(Vector{Point2f0}, collect(coordinates(rect)))
```
The element conversion is handled by `simplex_convert`, which also handles convert
between different face types:
```julia
decompose(QuadFace{Int}, rect) == convert(Vector{QuadFace{Int}}, collect(faces(rect)))
length(decompose(QuadFace{Int}, rect)) == 1
fs = decompose(GLTriangleFace, rect)
fs isa Vector{GLTriangleFace}
length(fs) == 2 # 2 triangles make up one quad ;)
```
`mesh` uses the most natural element type by default, which you can get with the unqualified Point type:
```julia
decompose(Point, rect) isa Vector{Point{2, Float64}}
```
You can also pass the element type to `mesh`:
```julia
m = GeometryBasics.mesh(rect, pointtype=Point2f0, facetype=QuadFace{Int})
```
You can also set the uv and normal type for the mesh constructor, which will then
calculate them for you, with the requested element type:
```julia
m = GeometryBasics.mesh(rect, uv=Vec2f0, normaltype=Vec3f0)
```

As you can see, the normals are automatically calculated,
the same is true for texture coordinates. You can overload this behavior by overloading
`normals` or `texturecoordinates` the same way as coordinates.
`decompose` works a bit different for normals/texturecoordinates, since they dont have their own element type.
Instead, you can use `decompose` like this:
```julia
decompose(UV(Vec2f0), rect)
decompose(Normal(Vec3f0), rect)
# the short form for the above:
decompose_uv(rect)
decompose_normals(rect)
```
You can also use `triangle_mesh`, `normal_mesh` and `uv_normal_mesh` to call the
`mesh` constructor with predefined element types (Point2/3f0, Vec2/3f0), and the requested attributes.
