# Point, Vec and Mat

GeometryBasics defines its own set of (small) Static Vectors and Matrices:
```julia
Point{N,T} <: StaticVector{N,T} <: AbstractVector{T}
Vec{N,T} <: StaticVector{N,T} <: AbstractVector{T}
Mat{Row, Column, T, L} <: AbstractMatrix{T}
```

These types are used throughout GeometryBasics to speed up calculations similar to how StaticArrays.jl does.

## Aliases

GeometryBasics exports common aliases for Point, Vec, Mat and Rect:

### Vec

|        |`T`(eltype) |`Float64` |`Float32` |`Int`     |`UInt`    |
|--------|------------|----------|----------|----------|----------|
|`N`(dim)|`Vec{N,T}`  |`Vecd{N}` |`Vecf{N}` |`Veci{N}` |`Vecui{N}`|
|`2`     |`Vec2{T}`   |`Vec2d`   |`Vec2f`   |`Vec2i`   |`Vec2ui`  |
|`3`     |`Vec3{T}`   |`Vec3d`   |`Vec3f`   |`Vec3i`   |`Vec3ui`  |

### Point

|        |`T`(eltype) |`Float64` |`Float32` |`Int`     |`UInt`    |
|--------|------------|----------|----------|----------|----------|
|`N`(dim)|`Point{N,T}`|`Pointd{N}`|`Pointf{N}`|`Pointi{N}`|`Pointui{N}`|
|`2`     |`Point2{T}` |`Point2d` |`Point2f` |`Point2i` |`Point2ui`|
|`3`     |`Point3{T}` |`Point3d` |`Point3f` |`Point3i` |`Point3ui`|

### Mat

|        |`T`(eltype) |`Float64` |`Float32` |`Int`     |`UInt`    |
|--------|------------|----------|----------|----------|----------|
|`N`(dim)|`Mat{N,T}`  |`Matd{N}` |`Matf{N}` |`Mati{N}` |`Matui{N}`|
|`2`     |`Mat2{T}`   |`Mat2d`   |`Mat2f`   |`Mat2i`   |`Mat2ui`  |
|`3`     |`Mat3{T}`   |`Mat3d`   |`Mat3f`   |`Mat3i`   |`Mat3ui`  |
