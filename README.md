
[![Build Status](https://api.travis-ci.org/JuliaGeometry/GeometryBasics.jl.svg?branch=master)](https://travis-ci.com/JuliaGeometry/GeometryBasics.jl)
[![Codecov](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl)

**Documentation**:
[![Docs - stable](https://img.shields.io/badge/docs-stable-lightgrey.svg)](http://juliageometry.github.io/GeometryBasics.jl/stable/)
[![Docs - dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliageometry.github.io/GeometryBasics.jl/dev)

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

Geometries can carry metadata:

```@repl quickstart
poi = meta(p1, city="Abuja", rainfall=1221.2)
```

Metadata is stored in a NamedTuple and can be retrieved as such:

```@repl quickstart
meta(poi)
```

Specific metadata attributes can be directly retrieved:

```@repl quickstart
poi.rainfall
```

To remove the metadata and keep only the geometry, use `metafree`:

```@repl quickstart
metafree(poi)
```

Geometries have predefined metatypes:

```@repl quickstart
multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)
```

Connect the points with lines:

```@repl quickstart
l1 = Line(p1, p2)
l2 = Line(p2, p3);
```

Connect the lines in a linestring:

```@repl quickstart
LineString([l1, l2])
```

Linestrings can also be constructed directly from points:

```@repl quickstart
LineString([p1, p2, p3])
```

The same goes for polygons:

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

Decompose the rectangle into four vertices:

```@repl quickstart
rect_vertices = decompose(Point{2, Float64}, rect)
```

Combine the vertices and faces into a triangle mesh:

```@repl quickstart
mesh = Mesh(rect_vertices, rect_faces)
```

Use `GeometryBasics.mesh` to get a mesh directly from a geometry:

```@repl quickstart
mesh = GeometryBasics.mesh(rect)
```
