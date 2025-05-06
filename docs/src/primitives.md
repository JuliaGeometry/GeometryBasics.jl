# Primitives

In GeometryBasics.jl, a `GeometryPrimitive` is an object from which a mesh can
be constructed.

## Existing GeometryPrimitives

GeometryBasics comes with a few predefined primitives:

#### HyperRectangle

A `Rect{D, T} = HyperRectangle{D, T}` is a D-dimensional axis-aligned
hyperrectangle defined by an origin and a size.

```@repl rects
using GeometryBasics
r1 = HyperRectangle{4, Float64}(Point{4, Float64}(0), Vec{4, Float64}(1))
r2 = Rect3f(Point3f(-1), Vec3f(2))
r3 = Rect2i(0, 0, 1, 1)
```

Rect2 supports normal and texture coordinate generation as well as tessellation.
Without tessellation, the coordinates of 2D Rects are defined in anti-clockwise order.
Rect3 supports normals and texture coordinates, but not tessellation.

Shorthands:

|        |`T`(eltype) |`Float64` |`Float32` |`Int`     |`UInt`    |
|--------|------------|----------|----------|----------|----------|
|`N`(dim)|`Rect{N,T}` |`Rectd{N}`|`Rectf{N}`|`Recti{N}`|`Rectui{N}`|
|`2`     |`Rect2{T}`  |`Rect2d`  |`Rect2f`  |`Rect2i`  |`Rect2ui` |
|`3`     |`Rect3{T}`  |`Rect3d`  |`Rect3f`  |`Rect3i`  |`Rect3ui` |

#### Sphere and Circle

`Circle` and `Sphere` are the 2 and 3 dimensional variants of `HyperSphere`.
They are defined by an origin and a radius.
While you can technically create a HyperSphere of any dimension, decomposition
is only defined in 2D and 3D.

```@setup hypersphere
using GeometryBasics
```
```@repl hypersphere
s1 = HyperSphere{4, Int}(Point{4, Int}(0), 5)
s2 = Sphere(Point3f(0, 0, 1), 1)
s3 = Circle(Point2d(0), 2.0)
```

Circle and Sphere support normal and texture coordinate generation as well as tessellation.
The coordinates of Circle are defined in anti-clockwise order.

#### Cylinder

A `Cylinder` is a 3D shape defined by two points and a radius.

```@setup cylinder
using GeometryBasics
```
```@repl cylinder
c = Cylinder(Point3f(-1, 0, 0), Point3f(0, 0, 1), 0.3f0) # start point, end point, radius
```

Cylinder supports normals and Tessellation, but currently no texture coordinates.

#### Cone

A `Cone` is also defined by two points and a radius, but the radius decreases to 0 from the start point to the tip.

```@setup cone
using GeometryBasics
```
```@repl cone
c = Cone(Point3f(-1, 0, 0), Point3f(0, 0, 1), 0.3f0) # start point, tip point, radius
```

Cone supports normals and Tessellation, but currently no texture coordinates.

#### Pyramid

`Pyramid` corresponds to a pyramid shape with a square base and four triangles
coming together into a sharp point.
It is defined by by the center point of the base, its height and its width.

```@setup pyramid
using GeometryBasics
```
```@repl pyramid
p = Pyramid(Point3f(0), 1f0, 0.3f0) # center, height, width
```

Pyramid supports normals, but currently no texture coordinates or tessellation

## Tessellation

In GeometryBasics `Tessellation` is a wrapper type for primitives which communicates
how dense the mesh generated from one should be.

```@setup tessellation
using GeometryBasics
```
```@repl tessellation
t = Tessellation(Cylinder(Point3f(0), Point3f(0,0,1), 0.2), 32) # 32 vertices for each circle
normal_mesh(t)

t = Tessellation(Rect2(Point2f(0), Vec2f(1)), (8, 6)) # 8 vertices in x direction by 6 in y direction
triangle_mesh(t)
```

## Primitive Interface / Implementing a new GeometryPrimitive

Every primitive should inherit from `GeometryPrimitive{Dim, eltype}` and implement at least `coordinates(primitive)` and `faces(primitive)` so that a mesh can be build from it.
This will also be enough to automatically generate normals for a 3D primitive and texture coordinates for a 2D primitive.
You can also implement functions to generate them directly with `normals(primitive)` and `texturecoordinates(primitive)`.
Depending on your primitive this might be necessary to get the normals and uvs you want.

To be compatible with `Tessellation` all of the functions mentioned above should implement a second tessellation argument.
This will be the second argument passed to the Tessellation constructor.
It's up to you to decide what makes sense here, though typically it's just an integer that more or less corresponds to the number of generated vertices.

#### Example

As an example, let's implement a parallelepiped, i.e. a 3D version or a parallelogram.
In this case we need an origin and 3 vectors telling us how far and in which directions the object extends.

```julia
struct Parallelepiped{T} <: GeometryPrimitive{3, T}
    origin::Point{3, T}
    v1::Vec{3, T}
    v2::Vec{3, T}
    v3::Vec{3, T}
end
```

Like the `Rect{3}`, this object comes with 8 unique positions which we want to return as its `coordinates`.

```julia
function GeometryBasics.coordinates(primitive::Parallelepiped{T}) where {T}
    o = primitive.origin
    v1 = primitive.v1; v2 = primitive.v2; v3 = primitive.v3
    return Point{3, T}[o, o+v2, o+v1+v2, o+v1,  o+v3, o+v2+v3, o+v1+v2+v3, o+v1+v3]
end
```

To connect these points into a mesh, we need to generate a set of faces.
The faces of a parallelepiped are parallelograms, which we can describe with `QuadFace`.
Here we should be conscious of the winding direction of faces.
They are often used to determine the front vs the backside of a (2D) face.
For example GeometryBasics normal generation and OpenGL's backface culling assume a counter-clockwise winding direction to correspond to a front-facing face.
This means that if we look at a face from outside the shape, the positions referred to by that face should be ordered counter-clockwise.
With that in mind the faces of our primitive become:

```julia
function GeometryBasics.faces(::Parallelepiped)
    return QuadFace{Int}[
        (1, 2, 3, 4), (5, 8, 7, 6), # facing -n3, +n3 (n3 being the normal of v1 x v2)
        (1, 5, 6, 2), (4, 3, 7, 8), # facing -n2, +n2
        (2, 6, 7, 3), (1, 4, 8, 5), # facing -n1, +n1
    ]
end
```

Note that you can check the correct winding direction fairly easily with Makie and the default generated normals.
After implementing faces and coordinates, you can create a mesh plot of your primitive with `Makie.mesh(primitive)`.
If the mesh reacts to light in a reasonable way, i.e. gets brighter when light shines on it, then your faces have the correct winding direction.
(It maybe useful to compare to other primitives like `Sphere(Point3f(0), 1f0)` here.)

Next on our TODO list are normals.
The default normals produced by `GeometryBasics.normal(primitive)` are vertex normals, which assume that a primitive to be smooth.
Since this is not the case for our primitive, we need to implement custom normals.
Here we could rely on `GeometryBasics.face_normal()` which returns a normal per face, but for this example we will implement them ourselves.

For our shape we want one normal per face, pointing in the normal direction of the corresponding 2D plane.
We can calculate the normal vector as `n = normalize(cross(v, w))` where v and w correspond to combinations of v1, v2 and v3.
To get them to act per face rather than per vertex, we need to overwrite the faces generated by `faces()`.
We can do that by creating a `FaceView` with a new set of faces which only act on normals.
Each of these new faces needs to refer to one normal by index to get what we want.

```julia
using LinearAlgebra
function GeometryBasics.normals(primitive::Parallelepiped)
    n1 = normalize(cross(primitive.v2, primitive.v3))
    n2 = normalize(cross(primitive.v3, primitive.v1))
    n3 = normalize(cross(primitive.v1, primitive.v2))
    ns = [-n3, n3, -n2, n2, -n1, n1]
    fs = QuadFace{Int}[1, 2, 3, 4, 5, 6] # = [QuadFace{Int}(1), QuadFace{Int}(2), ...]
    return FaceView(ns, fs)
end
```

As the last piece of the interface we can implement texture coordinates.
They generally refer to a 2D image with normalized 2D coordinates on a per-vertex basis.
There are many ways to define these coordinates.
Here we will partition the image in 2x3 even sized rectangular sections, split by the sign of the normal directions defined above.

```julia
function GeometryBasics.texturecoordinates(::Parallelepiped{T}) where {T}
    uvs = [Vec2f(x, y) for x in range(0, 1, length=4) for y in range(0, 1, 3)]
    fs = QuadFace{Int}[
        (1, 2, 5, 4),   (2, 3, 6, 5),
        (4, 5, 8, 7),   (5, 6, 9, 8),
        (7, 8, 11, 10), (8, 9, 12, 11)
    ]
    return FaceView(uvs, fs)
end
```
