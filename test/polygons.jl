@testset "Polygon" begin
    @testset "from points" begin
        points = connect([1, 2, 3, 4, 5, 6], Point2f)
        polygon = Polygon(points)
        @test polygon == Polygon(points)
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
