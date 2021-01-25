# Metadata

## Meta

The `Meta` method provides metadata handling capabilities in GeometryBasics.
Similarly to remove the metadata and keep only the geometry, use `metafree`, and
for vice versa i.e., remove the geometry and keep the metadata use `meta`.

### Syntax

```julia
meta(geometry, meta::NamedTuple)
meta(geometry; meta...)

metafree(meta-geometry)
meta(meta-geometry)
```

### Examples

```jldoctest meta
julia> using GeometryBasics

julia> p1 = Point(2.2, 3.6)
2-element Point{2,Float64} with indices SOneTo(2):
 2.2
 3.6

julia> poi = meta(p1, city="Abuja", rainfall=1221.2)
2-element PointMeta{2,Float64,Point{2,Float64},(:city, :rainfall),Tuple{String,Float64}} with indices SOneTo(2):
 2.2
 3.6
```

Metadata is stored in a NamedTuple and can be retrieved as such

```jldoctest meta
julia> meta(poi)
(city = "Abuja", rainfall = 1221.2)
```

Specific metadata attributes can be directly retrieved

```jldoctest meta
julia> poi.rainfall
1221.2

julia> metafree(poi)
2-element Point{2,Float64} with indices SOneTo(2):
 2.2
 3.6
```

For other geometries metatypes are predefined

```jldoctest meta
julia> multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)
1-element MultiPointMeta{Point{2,Float64},MultiPoint{2,Float64,Point{2,Float64},Array{Point{2,Float64},1}},(:city, :rainfall),Tuple{String,Float64}}:
 [2.2, 3.6]
```
In the above example we have also used geometry specific meta methods.


```jldoctest meta
julia> GeometryBasics.MetaType(Polygon)
PolygonMeta

julia> GeometryBasics.MetaType(Mesh)
MeshMeta
```
The metageometry objects are infact composed of the original geometry types.

```jldoctest meta
julia> GeometryBasics.MetaFree(PolygonMeta)
Polygon

julia> GeometryBasics.MetaFree(MeshMeta)
Mesh
```

## MetaT

In GeometryBasics we can a have tabular layout for a collection of meta-geometries
by putting them into a StructArray that extends the [Tables.jl](https://github.com/JuliaData/Tables.jl) API.

In practice it's not necessary for the geometry or metadata types to be consistent.
For example, a geojson format can have heterogeneous geometries. Hence, such cases require
automatic widening of the geometry data types to the most appropriate type.
The MetaT method works around the fact that, a collection of geometries and metadata
of different types can be represented tabularly whilst widening to the appropriate type.

### Syntax

```julia
MetaT(geometry, meta::NamedTuple)
MetaT(geometry; meta...)
```
Returns a `MetaT` that holds a geometry and its metadata `MetaT` acts the same as `Meta` method.
The difference lies in the fact that it is designed to handle geometries and metadata of different/heterogeneous types.

For example, while a Point MetaGeometry is a `PointMeta`, the MetaT representation is `MetaT{Point}`.

### Examples

```jldoctest meta
julia> MetaT(Point(1, 2), city = "Mumbai")
MetaT{Point{2,Int64},(:city,),Tuple{String}}([1, 2], (city = "Mumbai",))
```

For a tabular representation, an iterable of `MetaT` types can be passed on to a `metatable` method.

### Syntax

```julia
meta_table(iter)
```

### Examples

 Create an array of 2 linestrings:

```jldoctest meta
julia> ls = [LineString([Point(i, i+1), Point(i-1,i+5)]) for i in 1:2];

julia> coordinates.(ls)
2-element Array{Array{Point{2,Int64},1},1}:
 [[1, 2], [0, 6]]
 [[2, 3], [1, 7]]
```

Create a multi-linestring:

```jldoctest meta
julia> mls = MultiLineString(ls);

julia> coordinates.(mls)
2-element Array{Array{Point{2,Int64},1},1}:
 [[1, 2], [0, 6]]
 [[2, 3], [1, 7]]
```

Create a polygon:

```jldoctest meta
julia> poly = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)]);

julia> coordinates(poly)
4-element Array{Point{2,Int64},1}:
 [40, 40]
 [20, 45]
 [45, 30]
 [40, 40]
```

Put all geometries into an array:

```jldoctest meta
julia> geom = [ls..., mls, poly];
```

:Generate some random metadata:

```jldoctest meta
julia> prop = [(country_states = "India$(i)", rainfall = (i*9)/2) for i in 1:4]
4-element Array{NamedTuple{(:country_states, :rainfall),Tuple{String,Float64}},1}:
 (country_states = "India1", rainfall = 4.5)
 (country_states = "India2", rainfall = 9.0)
 (country_states = "India3", rainfall = 13.5)
 (country_states = "India4", rainfall = 18.0)

julia> feat = [MetaT(i, j) for (i,j) = zip(geom, prop)]; # create an array of MetaT
```

We can now generate a `StructArray` / `Table` with `meta_table`. Fields are accessed with `sa.main`, `sa.country_states`, `sa.rainfall`.

```jldoctest meta
julia> sa = meta_table(feat);
```

### Disadvantages

 * The MetaT is pretty generic in terms of geometry types, it's not subtype to
   geometries. eg : A `MetaT{Point, NamedTuple{Names, Types}}` is not subtyped to
   `AbstractPoint` like a `PointMeta` is.

 * This might cause problems on using `MetaT` with other constructors/methods
   inside or even outside GeometryBasics methods designed to work with the main `Meta` types.
