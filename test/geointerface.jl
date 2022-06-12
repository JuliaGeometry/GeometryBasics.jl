@testset "Basic types" begin
    point = Point(2, 3)
    GeoInterface.testgeometry(point)

    mp = MultiPoint([point, point])
    GeoInterface.testgeometry(mp)

    linestring = LineString(Point{2,Int}[(10, 10), (20, 20), (10, 40)])
    GeoInterface.testgeometry(linestring)

    multilinestring = MultiLineString([linestring, linestring])
    GeoInterface.testgeometry(multilinestring)

    poly = Polygon(rand(Point{2,Float32}, 5))
    GeoInterface.testgeometry(poly)

    triangle = Triangle(point, point, point)
    GeoInterface.testgeometry(triangle)

    polys = MultiPolygon([poly, poly])
    GeoInterface.testgeometry(polys)
end

@testset "Mesh" begin
    mesh = triangle_mesh(Sphere(Point3f(0), 1))
    GeoInterface.testgeometry(mesh)
end
