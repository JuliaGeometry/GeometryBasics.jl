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

        polygon = Polygon(points)
        mp = MultiPolygon([polygon, p1, p2])
        @test mp.polygons == [polygon, p1, p2]
        @test mp[1] == polygon
        @test mp[2] == p1
        @test size(mp) == (3,) # TODO: What does size even mean here?
        @test length(mp) == 3
        @test MultiPolygon(OffsetArray([polygon, p1, p2], 0)) == mp
        mp = MultiPolygon([Polygon(Point2f.()), p1, p2])

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

@testset "triangulation" begin
    tri = Triangle(Point2f(-0.5, -0.5), Point2f(0.8, 0), Point2f(0, 0.8))
    for phi in rand(0:2pi, 10)
        @test !in(Point2f(cos(phi), sin(phi)), tri)
    end
    for phi in rand(0:2pi, 10)
        @test in(0.2 * Point2f(cos(phi), sin(phi)), tri)
    end
    @test Point2f(0) in tri # sanity check

    # corner cases
    @test Point2f(-0.5) in tri
    @test Point2f(0, 0.8) in tri
    @test Point2f(0.8, 0) in tri

    # TODO: test snip directly?
    ps = Point2f[(0,0), (1,0), (1,1), (0,1)]
    @test decompose(GLTriangleFace, ps) == GLTriangleFace[(4,1,2), (2,3,4)]

    ps = [Point2f(cos(phi), sin(phi)) for phi in range(0, 2pi, length=8)[1:end-1]]
    fs = decompose(GLTriangleFace, ps)
    @test fs == GLTriangleFace[(7,1,2), (2,3,4), (4,5,6), (6,7,2), (2,4,6)]
    lfs = decompose(LineFace{Int32}, fs)
    for i in 1:7
        @test (LineFace(i, mod1(i+1, 7)) in lfs) || (LineFace(mod1(i+1, 7), i) in lfs)
    end

    @testset "earcut" begin
        ps = Point2i[(-1,-1), (1,-1), (1,1), (-1,1)]
        @test GeometryBasics.earcut_triangulate([Point{2, Int32}.(ps)]) == GLTriangleFace[(3,4,1), (1,2,3)]
        @test GeometryBasics.earcut_triangulate([Point{2, Int64}.(ps)]) == GLTriangleFace[(3,4,1), (1,2,3)]
        @test GeometryBasics.earcut_triangulate([Point{2, Float32}.(ps)]) == GLTriangleFace[(3,4,1), (1,2,3)]
        @test GeometryBasics.earcut_triangulate([Point{2, Float64}.(ps)]) == GLTriangleFace[(3,4,1), (1,2,3)]

        ps2 = Point2i[(0,-1), (1,0), (0,1), (-1,0)]
        @test faces(Polygon(Point2{Int32}.(ps), [Point2{Int32}.(ps2)])) == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
        @test faces(Polygon(Point2{Int64}.(ps), [Point2{Int64}.(ps2)])) == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
        @test faces(Polygon(Point2{Int8}.(ps),  [Point2{Int8}.(ps2)]))  == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
        @test faces(Polygon(Point2{Float32}.(ps), [Point2{Float32}.(ps2)])) == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
        @test faces(Polygon(Point2{Float64}.(ps), [Point2{Float64}.(ps2)])) == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
        @test faces(Polygon(Point2{Float16}.(ps), [Point2{Float16}.(ps2)])) == GLTriangleFace[(4,8,7), (5,8,1), (6,5,2), (3,7,6)]
    end
end