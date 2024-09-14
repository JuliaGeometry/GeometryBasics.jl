@testset "Polygon" begin
    @testset "from points" begin
        points = connect([1, 2, 3, 4, 5, 6], Point2f)
        polygon = Polygon(points)
        @test polygon == Polygon(points)
    end

    rect = Rect2f(0, 0, 1, 1)
    hole = Tesselation(Circle(Point2f(0.5), 0.2), 8)
    poly2 = Polygon(decompose(Point2f, rect), [decompose(Point2f, hole)])
    poly1 = Polygon(rect, [hole])
    @test poly1 == poly2
    @test poly1.exterior == decompose(Point2f, rect)
    @test poly1.interiors == [decompose(Point2f, hole)]

    fs = GLTriangleFace[
        (1, 9, 8), (5, 13, 12), (10, 9, 1), (4, 1, 8), (11, 10, 1), (3, 4, 8), 
        (12, 11, 1), (3, 8, 7), (12, 1, 2), (3, 7, 6), (5, 12, 2), (2, 3, 6), (6, 5, 2)]
    @test fs == faces(poly1)
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
end
