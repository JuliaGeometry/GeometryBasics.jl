# Meshes

```@setup 1
using Bonito
Page()
```

GeometryBasics defines two mesh types to work with - `Mesh` and `MetaMesh`

## Mesh

```@docs; canonical=false
Mesh
```

You can get data from a mesh using a few interface functions:
- `vertex_attributes(mesh) = mesh.vertex_attributes`
- `coordinates(mesh) = mesh.vertex_attributes[:position]`
- `normals(mesh) = mesh.vertex_attributes[:normal]`
- `texturecoordinates(mesh) = mesh.vertex_attributes[:uv]`
- `faces(mesh) = mesh.faces`

You can also grab the contents of `mesh.vertex_attributes` as if they were fields of the `Mesh`, e.g. `mesh.position` works.

### Custom Attributes

Meshes support arbitrary custom attributes beyond the standard position, normal, and UV coordinates. You can attach per-vertex or per-face data like material properties, identifiers, or computed values. These are stored in `mesh.vertex_attributes` and accessed using the same interface.

```julia
using GeometryBasics

# Define custom data type
struct Material
    emissivity::Float64
    absorptivity::Float64
    reflectivity::Float64
end

# Create mesh with custom attributes
points = [Point3f(0,0,0), Point3f(1,0,0), Point3f(0,1,0), Point3f(0,0,1)]
faces = [TriangleFace(1,2,3), TriangleFace(1,2,4), TriangleFace(1,3,4), TriangleFace(2,3,4)]

materials = [
    Material(0.1, 0.8, 0.1),
    Material(0.2, 0.7, 0.1),
    Material(0.0, 0.9, 0.1),
    Material(0.3, 0.6, 0.1)
]

face_names = ["bottom", "side1", "side2", "top"]

# Use per_face to create FaceViews for per-face attributes
mesh = GeometryBasics.mesh(
    points,
    faces,
    material=per_face(materials, faces),
    face_name=per_face(face_names, faces)
)

# Access custom attributes
mesh.material[2]  # Get material of second face
mesh.face_name[1]  # Get name of first face
```

This pattern is useful for physical simulations, rendering with material properties, or tagging mesh regions for analysis.

### FaceView


```@docs; canonical=false
FaceView
```

The purpose of FaceView is to allow you to add data that doesn't use the same vertex indices as `mesh.faces`
As a minimal example consider a mesh that is just one triangle, i.e. 3 position and one triangle face `TriangleFace(1,2,3)`.
Let's say we want to add a flat color to the triangle.
In this case we only have one color, but our face refers to 3 different vertices (3 different positions).
To avoid duplicating the color data, we can instead define a new triangle face `TriangleFace(1)` and add the color attribute as a `FaceView([color], [TriangleFace(1)])`.
If we ever need the mesh to be defined with just one common set of faces, i.e. no FaceView and appropriately duplicated vertex data, we can use `expand_faceviews(mesh)` to generate it.

On a larger scale this can be useful for memory and performance reason, e.g. when you do calculations with vertex attributes.
It can also simplify some definitions, like for example `Rect3`.
In that case we have 8 positions and 6 normals with FaceViews, or 24 without (assuming per-face normals).

For the relatively common case of per-face data, you can use the `per_face` convenience function.

```@docs
per_face
```

## MetaMesh

```@docs; canonical=false
MetaMesh
```

`MetaMesh` wraps a `Mesh` and allows you to attach global metadata that applies to the entire mesh rather than individual vertices or faces. This is useful for storing information like source file paths, transformation matrices, object identifiers, or simulation parameters.

```julia
using GeometryBasics

# Create a basic mesh
points = [Point3f(0,0,0), Point3f(1,0,0), Point3f(0,1,0)]
faces = [TriangleFace(1,2,3)]
mesh = GeometryBasics.mesh(points, faces; attribute=rand(3))

# Wrap with MetaMesh and add global metadata
meta_mesh = MetaMesh(mesh, source_file="model.obj", object_id=42, scale=1.5)

# Access metadata
meta_mesh[:source_file] # "model.obj"
meta_mesh[:object_id] # 42
meta_mesh.attribute # access vertex attributes via getproperty
# The underlying mesh is still accessible
meta_mesh.mesh
```

You can combine `MetaMesh` for global properties with per-face/per-vertex attributes for complete geometric and metadata representation:

```julia
# Create mesh with both per-face attributes and global metadata
mesh = GeometryBasics.mesh(
    points, faces,
    material=per_face(materials, faces),
    normal=face_normals(points, faces)
)

meta_mesh = MetaMesh(mesh, gltf_file="spacecraft.gltf", mass=150.0)
```

## How to create a mesh

### GeometryBasics

In GeometryBasics you mainly create meshes from primitives using a few constructors:
- `triangle_mesh(primitive)` generates the most basic mesh (i.e. positions and faces)
- `normal_mesh(primitive)` generates a mesh with normals (generated if the primitive doesn't implement `normal()`)
- `uv_mesh(primitive)` generates a mesh with texture coordinates (generated if the primitive doesn't implement `texturecoordinates()`)
- `uv_normal_mesh(primitive)` generates a mesh with normals and texture coordinates

Each of these constructors also includes keyword arguments for setting types, i.e. `pointtype`, `facetype`, `normaltype` and `uvtype` as appropriate.
Of course you can also construct a mesh directly from data, either with there various `Mesh()` or `GeometryBasics.mesh()` constructors.
The latter also include a `pointtype` and `facetype` conversion.

Finally there is also a `merge(::Vector{Mesh})` function which combines multiple meshes into a single one.
Note that this doesn't remove any data (e.g. hidden or duplicate vertices), and may remove `FaceView`s if they are incompatible between meshes.

### Meshing.jl

```@example 1
using Meshing
using GeometryBasics
using WGLMakie
using LinearAlgebra

gyroid(v) = cos(v[1])*sin(v[2])+cos(v[2])*sin(v[3])+cos(v[3])*sin(v[1])
gyroid_shell(v) = max(gyroid(v)-0.4,-gyroid(v)-0.4)
xr,yr,zr = ntuple(_->LinRange(0,pi*4,50),3)

A = [gyroid_shell((x,y,z)) for x in xr, y in yr, z in zr]
# generate directly using GeometryBasics API
# Rect specifies the sampling intervals
vts, fcs = isosurface(A, MarchingCubes())
# view with Makie
fcs = TriangleFace{Int}.(fcs)
vts = Point3d.(vts)
Makie.mesh(GeometryBasics.Mesh(vts, fcs), color=[norm(v) for v in vts])
```

### MeshIO.jl

The [`MeshIO.jl`](https://github.com/JuliaIO/MeshIO.jl) package provides load/save support for several file formats which store meshes.

```@example 1
using WGLMakie, GeometryBasics

m = Makie.loadasset("cat.obj")
Makie.mesh(m; color=Makie.loadasset("diffusemap.png"), axis=(; show_axis=false))
```
