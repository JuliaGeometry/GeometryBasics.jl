@testset "Polygon" begin
    @testset "from points" begin
        points = connect([1, 2, 3, 4, 5, 6], PointPoint2f(2))
        polygon = Polygon(points)
        @test polygon == Polygon(points)
    end

end

rect = Rect2f(0, 0, 1, 1)
hole = Circle(Point2f(0.5), 0.2)
poly2 = Polygon(decompose(Point2f, rect), [decompose(Point2f, hole)])
poly1 = Polygon(rect, [hole])
@test poly1 == poly2
@test poly.exterior == decompose(Point2f, rect)
@test poly.interiors == [decompose(Point2f, hole)]

faces(poly1)

GeometryBasics.earcut_triangulate([poly.exterior[[1, 2, 3, 4, 1]]])

poly = [
    [Point2f(100, 0), Point2f(100, 100), Point2f(0, 100), Point2f(0, 0)],
    # Following polylines define holes.
    [Point2f(75, 25), Point2f(75, 75), Point2f(25, 75), Point2f(25, 25)]
]

GeometryBasics.earcut_triangulate(poly)
