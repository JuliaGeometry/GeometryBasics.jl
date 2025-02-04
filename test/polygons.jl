@testset "Polygon" begin
    @testset "Constructors" begin
        points = connect([1, 2, 3, 4, 5, 6], Point2f)
        polygon = Polygon(points)
        @test polygon == Polygon(points)
        @test polygon == copy(polygon)
        @test coordinates(polygon) == points
        @test Polygon(OffsetArray(points)) == polygon

        interiors = [rand(Point2f, 4), rand(Point2f, 4)]
        exterior = rand(Point2f, 5)
        p1 = Polygon(exterior, interiors)
        @test p1.interiors == interiors
        @test p1.exterior == exterior
        p2 = Polygon(OffsetArray(exterior, 0), interiors)
        @test p2 == p1

        # TODO: promote polygon type automatically
        polygon = Polygon(Point2f.(points))
        mp = MultiPolygon([polygon, p1, p2])
        @test mp.polygons == [polygon, p1, p2]
        @test mp[1] == polygon
        @test mp[2] == p1
        @test size(mp) == (3,) # TODO: What does size even mean here?
        @test length(mp) == 3
        @test MultiPolygon(OffsetArray([polygon, p1, p2], 0)) == mp
    end

    rect = Rect2f(0, 0, 1, 1)
    hole = Tessellation(Circle(Point2f(0.5), 0.2), 8)
    poly2 = Polygon(decompose(Point2f, rect), [decompose(Point2f, hole)])
    poly1 = Polygon(rect, [hole])
    @test poly1 == poly2
    @test poly1.exterior == decompose(Point2f, rect)
    @test poly1.interiors == [decompose(Point2f, hole)]

    # triangulation is inconsistent...
    @test length(faces(poly1)) == 11
    ps = vcat(decompose(Point2f, rect), decompose(Point2f, hole))
    @test coordinates(poly1) == ps

    fs = GeometryBasics.earcut_triangulate([poly1.exterior[[1, 2, 3, 4, 1]]])
    @test fs == GLTriangleFace[(4,1,2), (2,3,4)]

    poly1 = [
        [Point2f(100, 0), Point2f(100, 100), Point2f(0, 100), Point2f(0, 0)],
        # Following polylines define holes.
        [Point2f(75, 25), Point2f(75, 75), Point2f(25, 75), Point2f(25, 25)]
    ]
    fs = GLTriangleFace[(4, 8, 7), (5, 8, 4), (3, 4, 7), (5, 4, 1), (2, 3, 7), (6, 5, 1), (2, 7, 6), (6, 1, 2)]
    @test fs == GeometryBasics.earcut_triangulate(poly1)

    # issue #236
    let p = SA[Point(0.0, 0.0), Point(1.0, 0.0), Point(0.0, 1.0)]
        @test Polygon(p) == Polygon(Array(p))
    end
end
