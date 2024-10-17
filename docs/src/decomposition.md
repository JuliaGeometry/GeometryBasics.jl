# Decomposition

## decompose functions

The `decompose` functions allow you to grab certain data from an `AbstractGeometry` like a mesh or primitive and convert it to a requested type, if possible.
They can also be used to convert an array of e.g. faces into a different face type directly.
The default decomposition implemented by GeoemtryBasics are:
- `decompose(::Type{<: Point}, source)` which collects data from `source` using `coordinates(source)` and converts it to the given point type.
- `decompose_normals([::Type{<: Vec},] source) = decompose([::Type{Normals{<: Vec}}},] source)` which collects data with `normals(source)` and converts it to the given Vec type.
- `decompose_uv([::Type{<: Vec},] source) = decompose([::Type{UV{<: Vec}}},] source)` which collects data with `texturecoordinates(source)` and converts it to the given Vec type. This function also exists with `UVW` texture coordinates.
- `decompose(::Type{<: AbstractFace}, source)` which collects data with `faces(source)` and converts it to the given face type. 

### Extending decompose

For `decompose` to work there needs to be a conversion from some element type to some target type. 
GeometryBasics relies on `GeometryBasics.convert_simplex(TargetType, value)` for this.
If you want to add new types to decompose, e.g. a new face type, you will need to add a method to that function.

## Primitive decomposition

GeometryBasics defines an interface to decompose geometry primitives into vertex attributes and faces.
The interface includes four functions:
- `coordinates(primitive[, nvertices])` which produces the positions associated with the primitive
- `faces(primitive[, nvertices])` which produces the faces which connect the vertex positions to a mesh
- `normals(primitive[, nvertices])` which optionally provide normal vectors of the primitive
- `texturecoordinates(primitive[, nvertices])` which optional provide texture coordinates (uv/uvw) of the primitive
