# Metadata
```julia
p1 = Point(2.2, 3.6)

# geometries can carry metadata
poi = meta(p1, city="Abuja", rainfall=1221.2)
2-element PointMeta{2,Int64,Point{2,Int64},(:city, :rainfall),Tuple{String,Float64}} with indices SOneTo(2):
 3
 1

# metadata is stored in a NamedTuple and can be retrieved as such
meta(poi)
(city = "Abuja", rainfall = 1221.2)

# specific metadata attributes can be directly retrieved
poi.rainfall
1221.2

# to remove the metadata and keep only the geometry, use metafree
metafree(poi)
2-element Point{2,Int64} with indices SOneTo(2):
 3
 1

# for other geometries metatypes are predefined
multipoi = MultiPointMeta([p1], city="Abuja", rainfall=1221.2)
1-element MultiPointMeta{Point{2,Int64},MultiPoint{2,Int64,Point{2,Int64},Array{Point{2,Int64},1}},(:city, :rainfall),Tuple{String,Float64}}:
[3, 1]
```

## MetaT
In GeometryBasics we can a have tabular layout for a collection meta-geometries by putting them into a StructArray that extends the Tables.jl API.

In practice it's not necessary for the geometry or metadata types to be consistent. Eg: A .geojson format can have heterogeneous geometries.
Hence, such cases require automatic widening of the geometry data types to the most appropriate type. The MetaT method works around the fact that, a collection of geometries and metadata of different types can be represented tabularly whilst widening to the appropriate type. 
### Syntax
```julia
    MetaT(geometry, meta::NamedTuple)
    MetaT(geometry; meta...)
```    
Returns a `MetaT` that holds a geometry and its metadata `MetaT` acts the same as `Meta` method.
The difference lies in the fact that it is designed to handle geometries and metadata of different/heterogeneous types.

eg: While a Point MetaGeometry is a `PointMeta`, the MetaT representation is `MetaT{Point}`

#### Example:
```julia
MetaT(Point(1, 2), city = "Mumbai")
MetaT{Point{2,Int64},(:city,),Tuple{String}}([1, 2], (city = "Mumbai",))
```

For a tabular representation, an iterable of `MetaT` types can be passed on to a `metatable` method.

### Syntax
```julia
    meta_table(iter)
```    
#### Example:
```julia
using DataFrames
# Create an array of 2 linestrings 
ls = [LineString([Point(i, i+1), Point(i-1,i+5)]) for i in 1:2]

# Create a MultiLineString 
mls = MultiLineString(ls)

# Create a Polygon
poly = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)])

# Put all of it in an Array
geom = [ls..., mls, poly]
    
# Generate some random metadata
prop = [(country_states = "India$(i)", rainfall = (i*9)/2) for i in 1:4]
    
# Create an Array of MetaT
feat = [MetaT(i, j) for (i,j) = zip(geom, prop)]

# Generate a StructArray/Table
sa = meta_table(feat)

sa.main
sa.country_states
sa.rainfall
```

### Disadvantages:
 * The MetaT is pretty generic in terms of geometry types, it's not subtype to geometries. 
 eg : A `MetaT{Point, NamedTuple{Names, Types}}` is not subtyped to `AbstractPoint` like a `PointMeta` is.
 * This might cause problems on using `MetaT` with other constructors/methods inside or even outside GeometryBasics methods designed to work with the main `Meta` types.

