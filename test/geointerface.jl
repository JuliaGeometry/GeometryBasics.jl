@testset "Basic types" begin
    point = Point(2, 3)
    GeoInterface.testgeometry(point)
    @test GeoInterface.ncoord(point) == 2
    @test GeoInterface.getcoord(point, 2) == 3
    @test GeoInterface.coordinates(point) == [2, 3]

    mp = MultiPoint([point, point])
    GeoInterface.testgeometry(mp)
    @test GeoInterface.ngeom(mp) == 2
    @test GeoInterface.getgeom(mp, 2) == point
    @test GeoInterface.coordinates(mp) == [[2, 3], [2, 3]]

    linestring = LineString(Point{2,Int}[(10, 10), (20, 20), (10, 40)])
    GeoInterface.testgeometry(linestring)
    @test GeoInterface.ngeom(linestring) == 3
    @test GeoInterface.getgeom(linestring, 1) == Point(10, 10)
    @test GeoInterface.getgeom(linestring, 2) == Point(20, 20)
    @test GeoInterface.getgeom(linestring, 3) == Point(10, 40)
    @test GeoInterface.coordinates(linestring) == [[10, 10], [20, 20], [10, 40]]

    multilinestring = MultiLineString([linestring, linestring])
    GeoInterface.testgeometry(multilinestring)
    @test GeoInterface.coordinates(multilinestring) ==
          [[[10, 10], [20, 20], [10, 40]], [[10, 10], [20, 20], [10, 40]]]

    poly = Polygon(rand(Point{2,Float32}, 5), [rand(Point{2,Float32}, 5)])
    GeoInterface.testgeometry(poly)
    @test length(GeoInterface.coordinates(poly)) == 2
    @test length(GeoInterface.coordinates(poly)[1]) == 5

    triangle = Triangle(point, point, point)
    GeoInterface.testgeometry(triangle)
    @test length(GeoInterface.coordinates(triangle)) == 1
    @test length(GeoInterface.coordinates(triangle)[1]) == 3

    polys = MultiPolygon([poly, poly])
    GeoInterface.testgeometry(polys)
    @test length(GeoInterface.coordinates(polys)) == 2
    @test length(GeoInterface.coordinates(polys)[1]) == 2
    @test length(GeoInterface.coordinates(polys)[1][1]) == 5
end

@testset "Mesh" begin
    mesh = triangle_mesh(Sphere(Point3f(0), 1))
    GeoInterface.testgeometry(mesh)
end
