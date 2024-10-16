# Meshes

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

### FaceView


```@docs; canoniical=false
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


## MetaMesh

```julia; canonical=false
MetaMesh
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

### MeshIO.jl

The [`MeshIO.jl`](https://github.com/JuliaIO/MeshIO.jl) package provides load/save support for several file formats which store meshes.

```@example
using GLMakie, GLMakie.FileIO, GeometryBasics

m = load(GLMakie.assetpath("cat.obj"))
GLMakie.mesh(m; color=load(GLMakie.assetpath("diffusemap.png")), axis=(; show_axis=false))
```
