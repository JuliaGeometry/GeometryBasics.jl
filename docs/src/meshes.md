# Meshes

A mesh consists of a sequence of [`Polytope`](@ref)s. For example, in a 3D dimensional space, a mesh would be a sequence of triangles (2-simplexes).

## Types

* [`AbstractMesh`](@ref)
* [`Mesh`](@ref)

## How to create a mesh

To create a mesh one can provide one the following.
* A list of points and faces.
* A list of polytopes.

First, let's create four points and four faces. Each face is an integer connecting the points according to their array index.

```jldoctest
julia> mypoints = [
           Point3f(0,0,0),
           Point3f(0,0,1),
           Point3f(0,1,0),
           Point3f(1,0,0)
       ]
4-element Vector{Point{3, Float32}}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 1.0]
 [0.0, 1.0, 0.0]
 [1.0, 0.0, 0.0]

julia> myfaces = [
           TriangleFace(1,2,3),
           TriangleFace(1,2,4),
           TriangleFace(1,3,4),
           TriangleFace(2,3,4)
       ]
4-element Vector{TriangleFace{Int64}}:
 TriangleFace(1, 2, 3)
 TriangleFace(1, 2, 4)
 TriangleFace(1, 3, 4)
 TriangleFace(2, 3, 4)


julia> mymesh = Mesh(mypoints, myfaces)
Mesh{3, Float32, Triangle}:
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])
```

As seen above, the mesh is just a sequence of triangles. Next, let's create a similar `Mesh` by providing the triangles directly.

```jldoctest
julia> mytriangles = [
           Triangle(pts[[1,2,3]]...),
           Triangle(pts[[1,2,4]]...),
           Triangle(pts[[1,3,4]]...),
           Triangle(pts[[2,3,4]]...)
       ]
4-element Vector{GeometryBasics.Ngon{3, Float32, 3, Point{3, Float32}}}:
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])

julia> mymesh2 = Mesh(mytriangles)
Mesh{3, Float32, Triangle}:
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 0.0, 1.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 0.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])
 Triangle(Float32[0.0, 0.0, 1.0], Float32[0.0, 1.0, 0.0], Float32[1.0, 0.0, 0.0])
```


### Meshing.jl

### MeshIO.jl

The [`MeshIO.jl`](https://github.com/JuliaIO/MeshIO.jl) package provides load/save support for several file formats which store meshes.

## How to access data

The following functions can be called on an [`AbstractMesh`](@ref) to access its underlying data.

* [`faces`](@ref)
* [`coordinates`](@ref)
* [`texturecoordinates`](@ref)
* [`normals`](@ref)


```jldoctest
julia> GeometryBasics.faces(mymesh)
4-element Vector{TriangleFace{Int64}}:
 TriangleFace(1, 2, 3)
 TriangleFace(1, 2, 4)
 TriangleFace(1, 3, 4)
 TriangleFace(2, 3, 4)

julia> GeometryBasics.coordinates(mymesh)
4-element Vector{Point{3, Float32}}:
 [0.0, 0.0, 0.0]
 [0.0, 0.0, 1.0]
 [0.0, 1.0, 0.0]
 [1.0, 0.0, 0.0]
```

Note that these functions may not apply to all meshes. For example, `mymesh2`
above was not created with `Faces` so `faces` will return `nothing`.

```jldoctest
julia> GeometryBasics.faces(mymesh2)

```