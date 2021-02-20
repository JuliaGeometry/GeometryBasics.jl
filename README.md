
[![Build Status](https://api.travis-ci.org/JuliaGeometry/GeometryBasics.jl.svg?branch=master)](https://travis-ci.com/JuliaGeometry/GeometryBasics.jl)
[![Codecov](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl)

**Documentation**:
[![Docs - stable](https://img.shields.io/badge/docs-stable-lightgrey.svg)](http://juliageometry.github.io/GeometryBasics.jl/stable/)
[![Docs - dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliageometry.github.io/GeometryBasics.jl/dev)

# GeometryBasics.jl

Basic Geometry Types.
This package aims to offer a standard set of Geometry types, which easily work with metadata, query frameworks on geometries and different memory layouts.
The aim is to create a solid basis for Graphics/Plotting, finite elements analysis, Geo applications, and general geometry manipulations - while offering a julian API, that still allows performant C-interop.

This package is a replacement for the discontinued [GeometryTypes](https://github.com/JuliaGeometry/GeometryTypes.jl/).

## Quick start

```julia
using GeometryBasics

# create some points
julia> p1 = Point(3, 1)
2-element Point{2,Int64} with indices SOneTo(2):
 3
 1

julia> p2 = Point(1, 3);

julia> p3 = Point(4, 4);

# geometries can carry metadata
julia> poi = meta(p1, city="Abuja", rainfall=1221.2)
2-element PointMeta{2,Int64,Point{2,Int64},(:city, :rainfall),Tuple{String,Float64}} with indices SOneTo(2):
 3
 1

# metadata is stored in a NamedTuple and can be retrieved as such
julia> meta(poi)
(city = "Abuja", rainfall = 1221.2)

# specific metadata attributes can be directly retrieved
julia> poi.rainfall
1221.2

# to remove the metadata and keep only the geometry, use metafree
julia> metafree(poi)
2-element Point{2,Int64} with indices SOneTo(2):
 3
 1

# for other geometries metatypes are predefined
julia> multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)
1-element MultiPointMeta{Point{2,Int64},MultiPoint{2,Int64,Point{2,Int64},Array{Point{2,Int64},1}},(:city, :rainfall),Tuple{String,Float64}}:
[3, 1]

# connect the points with lines
julia> l1 = Line(p1, p2)
Line([3, 1] => [1, 3])

julia> l2 = Line(p2, p3);

# connect the lines in a linestring
julia> LineString([l1, l2])
2-element LineString{2,Int64,Point{2,Int64},Array{GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}},1}}:
 Line([3, 1] => [1, 3])
 Line([1, 3] => [4, 4])

# linestrings can also be constructed directly from points
julia> LineString([p1, p2, p3])
2-element LineString{2,Int64,Point{2,Int64},Base.ReinterpretArray{GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}},1,Tuple{Point{2,Int64},Point{2,Int64}},TupleView{Tuple{Point{2,Int64},Point{2,Int64}}, 1}}}:
 Line([3, 1] => [1, 3])
 Line([1, 3] => [4, 4])

# the same goes for polygons
julia> Polygon(Point{2, Int}[(3, 1), (4, 4), (2, 4), (1, 2), (3, 1)])
Polygon{2,Int64,Point{2,Int64},LineString{2,Int64,Point{2,Int64},Base.ReinterpretArray{GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}},1,Tuple{Point{2,Int64},Point{2,Int64}},TupleView{Tuple{Point{2,Int64},Point{2,Int64}}, 1}}},Array{LineString{2,Int64,Point{2,Int64},Base.ReinterpretArray{GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}},1,Tuple{Point{2,Int64},Point{2,Int64}},TupleView{Tuple{Point{2,Int64},Point{2,Int64}}, 1}}},1}}(GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}}[Line([3, 1] => [4, 4]), Line([4, 4] => [2, 4]), Line([2, 4] => [1, 2]), Line([1, 2] => [3, 1])], LineString{2,Int64,Point{2,Int64},Base.ReinterpretArray{GeometryBasics.Ngon{2,Int64,2,Point{2,Int64}},1,Tuple{Point{2,Int64},Point{2,Int64}},TupleView{Tuple{Point{2,Int64},Point{2,Int64}}, 1}}}[])

# create a rectangle placed at the origin with unit widths
julia> rect = Rect(Vec(0.0, 0.0), Vec(1.0, 1.0))
GeometryBasics.HyperRectangle{2,Float64}([0.0, 0.0], [1.0, 1.0])

# decompose the rectangle into two triangular faces
julia> rect_faces = decompose(TriangleFace{Int}, rect)
2-element Array{NgonFace{3,Int64},1}:
 TriangleFace(1, 2, 4)
 TriangleFace(1, 4, 3)

# decompose the rectangle into four vertices
julia> rect_vertices = decompose(Point{2, Float64}, rect)
4-element Array{Point{2,Float64},1}:
 [0.0, 0.0]
 [1.0, 0.0]
 [0.0, 1.0]
 [1.0, 1.0]

# combine the vertices and faces into a triangle mesh
julia> mesh = Mesh(rect_vertices, rect_faces)
Mesh{2, Float64, Triangle}:
 Triangle([0.0, 0.0], [1.0, 0.0], [1.0, 1.0])
 Triangle([0.0, 0.0], [1.0, 1.0], [0.0, 1.0])
```
