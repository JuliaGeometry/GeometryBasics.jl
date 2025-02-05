# Bounding Boxes

You can generate an axis aligned bounding box for any `AbstractGeometry` by calling `Rect(geom)`.
Depending on the object this will either rely on `coordinates(geom)` or a specialized method.
You can also create a bounding box of set dimension or type by adding the related parameters.


```@repl
using GeometryBasics

s = Circle(Point2f(0), 1f0)
Rect(s) # specialized, exact bounding box
Rect3(s)
Rect3d(s)
RectT{Float64}(s)
Rect(GeometryBasics.mesh(s)) # using generated coordinates in mesh
```

## Extending

If you want to add a specialized bounding box method you should implement `Rect{N, T}(geom) = ...`.
All other methods funnel into that one, defaulting to the same `N, T` that the given `AbstractGeometry{N, T}` has.
GeometryBasics allows the user given dimension `N` to be smaller or equal to that of the geometry.
This is checked with `GeometryBasics.bbox_dim_check(user_dim, geom_dim)` which you may reuse.

```@example
function Rect{N, T}(a::HyperSphere{N2}) where {N, N2, T}
    bbox_dim_check(N, N2)
    return Rect{N, T}(minimum(a), widths(a))
end
```