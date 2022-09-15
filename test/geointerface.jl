@testset "Basic types" begin
    point = Point(2, 3)
    @test testgeometry(point)
    @test ncoord(point) == 2
    @test getcoord(point, 2) == 3
    @test GeoInterface.coordinates(point) == [2, 3]

    mp = MultiPoint([point, point])
    @test testgeometry(mp)
    @test ngeom(mp) == 2
    @test getgeom(mp, 2) == point
    @test GeoInterface.coordinates(mp) == [[2, 3], [2, 3]]

    linestring = LineString(Point{2,Int}[(10, 10), (20, 20), (10, 40)])
    @test testgeometry(linestring)
    @test ngeom(linestring) == 3
    @test getgeom(linestring, 1) == Point(10, 10)
    @test getgeom(linestring, 2) == Point(20, 20)
    @test getgeom(linestring, 3) == Point(10, 40)
    @test GeoInterface.coordinates(linestring) == [[10, 10], [20, 20], [10, 40]]

    multilinestring = MultiLineString([linestring, linestring])
    @test testgeometry(multilinestring)
    @test GeoInterface.coordinates(multilinestring) ==
          [[[10, 10], [20, 20], [10, 40]], [[10, 10], [20, 20], [10, 40]]]

    poly = Polygon(rand(Point{2,Float32}, 5), [rand(Point{2,Float32}, 5)])
    @test testgeometry(poly)
    @test length(GeoInterface.coordinates(poly)) == 2
    @test length(GeoInterface.coordinates(poly)[1]) == 5

    triangle = Triangle(point, point, point)
    @test testgeometry(triangle)
    @test length(GeoInterface.coordinates(triangle)) == 1
    @test length(GeoInterface.coordinates(triangle)[1]) == 3

    polys = MultiPolygon([poly, poly])
    @test testgeometry(polys)
    @test length(GeoInterface.coordinates(polys)) == 2
    @test length(GeoInterface.coordinates(polys)[1]) == 2
    @test length(GeoInterface.coordinates(polys)[1][1]) == 5
end

@testset "Mesh" begin
    mesh = triangle_mesh(Sphere(Point3f(0), 1))
    @test testgeometry(mesh)
end

@testset "Convert" begin
    # convert GeoJSON geometry types to GeometryBasics via the GeoInterface
    point_str = """{"type":"Point","coordinates":[30.1,10.1]}"""
    point_3d_str = """{"type":"Point","coordinates":[30.1,10.1,5.1]}"""
    linestring_str = """{"type":"LineString","coordinates":[[30.1,10.1],[10.1,30.1],[40.1,40.1]]}"""
    polygon_str = """{"type":"Polygon","coordinates":[[[30.1,10.1],[40.1,40.1],[20.1,40.1],[10.1,20.1],[30.1,10.1]]]}"""
    polygon_hole_str = """{"type":"Polygon","coordinates":[[[35.1,10.1],[45.1,45.1],[15.1,40.1],[10.1,20.1],[35.1,10.1]],[[20.1,30.1],[35.1,35.1],[30.1,20.1],[20.1,30.1]]]}"""
    multipoint_str = """{"type":"MultiPoint","coordinates":[[10.1,40.1],[40.1,30.1],[20.1,20.1],[30.1,10.1]]}"""
    multilinestring_str = """{"type":"MultiLineString","coordinates":[[[10.1,10.1],[20.1,20.1],[10.1,40.1]],[[40.1,40.1],[30.1,30.1],[40.1,20.1],[30.1,10.1]]]}"""
    multipolygon_str = """{"type":"MultiPolygon","coordinates":[[[[30.1,20.1],[45.1,40.1],[10.1,40.1],[30.1,20.1]]],[[[15.1,5.1],[40.1,10.1],[10.1,20.1],[5.1,10.1],[15.1,5.1]]]]}"""
    multipolygon_hole_str = """{"type":"MultiPolygon","coordinates":[[[[40.1,40.1],[20.1,45.1],[45.1,30.1],[40.1,40.1]]],[[[20.1,35.1],[10.1,30.1],[10.1,10.1],[30.1,5.1],[45.1,20.1],[20.1,35.1]],[[30.1,20.1],[20.1,15.1],[20.1,25.1],[30.1,20.1]]]]}"""

    point_json = GeoJSON.read(point_str)
    point_3d_json = GeoJSON.read(point_3d_str)
    linestring_json = GeoJSON.read(linestring_str)
    polygon_json = GeoJSON.read(polygon_str)
    polygon_hole_json = GeoJSON.read(polygon_hole_str)
    multipoint_json = GeoJSON.read(multipoint_str)
    multilinestring_json = GeoJSON.read(multilinestring_str)
    multipolygon_json = GeoJSON.read(multipolygon_str)
    multipolygon_hole_json = GeoJSON.read(multipolygon_hole_str)

    point_gb = GeoInterface.convert(Point, point_json)
    point_3d_gb = GeoInterface.convert(Point, point_3d_json)
    linestring_gb = GeoInterface.convert(LineString, linestring_json)
    polygon_gb = GeoInterface.convert(Polygon, polygon_json)
    polygon_hole_gb = GeoInterface.convert(Polygon, polygon_hole_json)
    multipoint_gb = GeoInterface.convert(MultiPoint, multipoint_json)
    multilinestring_gb = GeoInterface.convert(MultiLineString, multilinestring_json)
    multipolygon_gb = GeoInterface.convert(MultiPolygon, multipolygon_json)
    multipolygon_hole_gb = GeoInterface.convert(MultiPolygon, multipolygon_hole_json)

    @test point_gb === Point{2, Float64}(30.1, 10.1)
    @test point_3d_gb === Point{3, Float64}(30.1, 10.1, 5.1)
    @test linestring_gb isa LineString
    # TODO, what should we do exactly with linestrings?
    # @test length(linestring_gb) == 2
    # @test eltype(linestring_gb) == Line{2, Float64}
    @test polygon_gb isa Polygon
    @test isempty(polygon_gb.interiors)
    @test polygon_hole_gb isa Polygon
    @test length(polygon_hole_gb.interiors) == 1
    @test multipoint_gb isa MultiPoint
    @test length(multipoint_gb) == 4
    @test multipoint_gb[4] === Point{2, Float64}(30.1, 10.1)
    @test multilinestring_gb isa MultiLineString
    @test length(multilinestring_gb) == 2
    @test multipolygon_gb isa MultiPolygon
    @test length(multipolygon_gb) == 2
    @test multipolygon_hole_gb isa MultiPolygon
    @test length(multipolygon_hole_gb) == 2
    @test length(multipolygon_hole_gb[1].interiors) == 0
    @test length(multipolygon_hole_gb[2].interiors) == 1
end
