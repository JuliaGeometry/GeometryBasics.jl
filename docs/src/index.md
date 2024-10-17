# GeometryBasics.jl

Basic geometry types.

This package aims to offer a standard set of geometry types that easily work
with metadata, query frameworks on geometries and different memory layouts. The
aim is to create a solid basis for graphics/plotting, finite element analysis,
geo applications, and general geometry manipulations - while offering a Julian
API that still allows performant C-interop.

This package is a replacement for the discontinued [GeometryTypes](https://github.com/JuliaGeometry/GeometryTypes.jl/).

## Quick start

Create some points:

```@repl quickstart
using GeometryBasics

p1 = Point(3, 1)
p2 = Point(1, 3);
p3 = Point(4, 4);
```

Connect pairs of points as line segments:

```@repl quickstart
l1 = Line(p1, p2)
l2 = Line(p2, p3);
```

Or connect multiple points as a linestring:

```@repl quickstart
LineString([p1, p2, p3])
```

You can also create polygons from points:

```@repl quickstart
Polygon(Point{2, Int}[(3, 1), (4, 4), (2, 4), (1, 2), (3, 1)])
```

Create a rectangle placed at the origin with unit width and height:

```@repl quickstart
rect = Rect(Vec(0.0, 0.0), Vec(1.0, 1.0))
```

Decompose the rectangle into two triangular faces:

```@repl quickstart
rect_faces = decompose(TriangleFace{Int}, rect)
```

Decompose the rectangle into four positions:

```@repl quickstart
rect_positions = decompose(Point{2, Float64}, rect)
```

Combine the vertices and faces into a triangle mesh:

```@repl quickstart
mesh = Mesh(rect_positions, rect_faces)
```

Use `GeometryBasics.mesh` to get a mesh directly from a geometry:

```@repl quickstart
mesh = GeometryBasics.mesh(rect)
```
