# Decomposition


## Displaying primitives

To display geometry primitives, they need to be decomposable.
This can be done for any arbitrary primitive, by overloading the following interface:

```julia
# Let's take SimpleRectangle as an example:
# Below is a minimal set of decomposable attributes to build up a triangle mesh:
isdecomposable(::Type{T}, ::Type{HR}) where {T<:Point, HR<:SimpleRectangle} = true
isdecomposable(::Type{T}, ::Type{HR}) where {T<:Face, HR<:SimpleRectangle} = true

# This is an example implementation of `decompose` for points.
function GeometryBasics.decompose(P::Type{Point{3, PT}}, r::SimpleRectangle, resolution=(2,2)) where PT
    w,h = resolution
    vec(
        PT[
            (x,y,0)
            for x in range(r.x, stop = r.x+r.w, length = w),
                y in range(r.y, stop = r.y+ r .h, length = h)
        ]
    )
end

function GeometryBasics.decompose(::Type{T}, r::SimpleRectangle, resolution=(2,2)) where T <: Face
    w,h = resolution
    Idx = LinearIndices(resolution)
    faces = vec([Face{4, Int}(
            Idx[i, j], Idx[i+1, j],
            Idx[i+1, j+1], Idx[i, j+1]
        ) for i=1:(w-1), j=1:(h-1)]
    )
    decompose(T, faces)
end
```

With these methods defined, this constructor will magically work:

```julia
rect = SimpleRectangle(0, 0, 1, 1)
m = GLNormalMesh(rect)
vertices(m) == decompose(Point3f0, rect)

faces(m) == decompose(GLTriangle, rect) # GLFace{3} == GLTriangle
normals(m) # automatically calculated from mesh
```

As you can see, the normals are automatically calculated only with the faces and points.
You can overwrite that behavior by also defining decompose for the `Normal` type!
