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

```@repl meta
using GeometryBasics
p1 = Point(2.2, 3.6)
poi = meta(p1, city="Abuja", rainfall=1221.2)
```

Metadata is stored in a NamedTuple and can be retrieved as such:

```@repl meta
meta(poi)
```

Specific metadata attributes can be directly retrieved:

```@repl meta
poi.rainfall
metafree(poi)
```

Metatypes are predefined for geometries:

```@repl meta
multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)
```

(In the above example we have also used a geometry-specific meta method.)

```@repl meta
GeometryBasics.MetaType(Polygon)
GeometryBasics.MetaType(Mesh)
```

The metageometry objects are infact composed of the original geometry types.

```@repl meta
GeometryBasics.MetaFree(PolygonMeta)
GeometryBasics.MetaFree(MeshMeta)
```

## MetaT

In GeometryBasics we can have tabular layout for a collection of meta-geometries
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

```@repl meta
MetaT(Point(1, 2), city = "Mumbai")
```

For a tabular representation, an iterable of `MetaT` types can be passed on to a `meta_table` method.

### Syntax

```julia
meta_table(iter)
```

### Examples

 Create an array of 2 linestrings:

```@repl meta
ls = [LineString([Point(i, i+1), Point(i-1,i+5)]) for i in 1:2];
coordinates.(ls)
```

Create a multi-linestring:

```@repl meta
mls = MultiLineString(ls);
coordinates.(mls)
```

Create a polygon:

```@repl meta
poly = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)]);
coordinates(poly)
```

Put all geometries into an array:

```@repl meta
geom = [ls..., mls, poly];
```

Generate some random metadata:

```@repl meta
prop = [(country_states = "India$(i)", rainfall = (i*9)/2) for i in 1:4]
feat = [MetaT(i, j) for (i,j) = zip(geom, prop)]; # create an array of MetaT
```

We can now generate a `StructArray` / `Table` with `meta_table`:

```@repl meta
sa = meta_table(feat);
```

The data can be accessed through `sa.main` and the metadata through
`sa.country_states` and `sa.rainfall`. Here we print only the type names of the
data items for brevity:

```@repl meta
[nameof.(typeof.(sa.main)) sa.country_states sa.rainfall]
```

### Disadvantages

 * The MetaT is pretty generic in terms of geometry types, it's not subtype to
   geometries. eg : A `MetaT{Point, NamedTuple{Names, Types}}` is not subtyped to
   `AbstractPoint` like a `PointMeta` is.

 * This might cause problems on using `MetaT` with other constructors/methods
   inside or even outside GeometryBasics methods designed to work with the main `Meta` types.
