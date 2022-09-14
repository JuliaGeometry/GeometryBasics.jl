using GeometryBasicsCore
using LinearAlgebra
using Test

@testset "Point & Vec type" begin
    @testset "conversion" begin
        @test convert(Point, (2, 3)) === Point(2, 3)
        @test convert(Point, (2.0, 3)) === Point(2.0, 3.0)
    end

    @testset "broadcast" begin
        @testset for T in (Vec, Point)
            x = [T(2, 3), T(7, 3)]

            @test [T(4, 9), T(14, 9)] == x .* T(2, 3)
            @test [T(4, 6), T(9, 6)] == x .+ T(2, 3)
            @test [T(0, 0), T(5, 0)] == x .- T(2, 3)
        end
    end
end

@testset "constructors" begin
    @testset "LineFace" begin
        points = connect([1, 2, 3, 4, 5, 6], Point{2})
        linestring = LineString(points)
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        points = rand(Point{2, Float64}, 4)
        linestring = LineString(points, 2)
        @test linestring == [Line(points[1], points[2]), Line(points[3], points[4])]

        linestring = LineString([points[1] => points[2], points[2] => points[3]])
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        faces = [1, 2, 3]
        linestring = LineString(points, faces)
        @test linestring == LineString([points[1] => points[2], points[2] => points[3]])
        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        linestring = LineString(points, faces, 2)
        @test linestring == LineString([a => b, c => d])

        faces = [LineFace(1, 2)
        , LineFace(3, 4)]
        linestring = LineString(points, faces)
        @test linestring == LineString([a => b, c => d])
    end

    @testset "Polygon" begin
        points = connect([1, 2, 3, 4, 5, 6], Point{2})
        polygon = Polygon(points)
        @test polygon == Polygon(LineString(points))

        points = rand(Point{2, Float64}, 4)
        linestring = LineString(points, 2)
        @test Polygon(points, 2) == Polygon(linestring)

        faces = [1, 2, 3]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))

        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        polygon = Polygon(points, faces, 2)
        @test polygon == Polygon(LineString(points, faces, 2))

        faces = [LineFace(1, 2), LineFace(3, 4)]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))
        @test ndims(polygon) === 2
    end

    @testset "Multi geometries" begin
        # coordinates from https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry#Geometric_objects
        points = Point{2, Int}[(10, 40), (40, 30), (20, 20), (30, 10)]
        multipoint = MultiPoint(points)
        @test size(multipoint) === size(points)
        @test multipoint[3] === points[3]

        linestring1 = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2, Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        @test size(multilinestring) === (2,)
        @test multilinestring[1] === linestring1
        @test multilinestring[2] === linestring2

        polygon11 = Polygon(Point{2, Int}[(30, 20), (45, 40), (10, 40), (30, 20)])
        polygon12 = Polygon(Point{2, Int}[(15, 5), (40, 10), (10, 20), (5, 10), (15, 5)])
        multipolygon1 = MultiPolygon([polygon11, polygon12])
        @test size(multipolygon1) === (2,)
        @test multipolygon1[1] === polygon11
        @test multipolygon1[2] === polygon12

        polygon21 = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)])
        polygon22 = Polygon(LineString(Point{2, Int}[(20, 35), (10, 30), (10, 10), (30, 5), (45, 20), (20, 35)]),
            [LineString(Point{2, Int}[(30, 20), (20, 15), (20, 25), (30, 20)])])
        multipolygon2 = MultiPolygon([polygon21, polygon22])
        @test size(multipolygon2) === (2,)
        @test multipolygon2[1] === polygon21
        @test multipolygon2[2] === polygon22
    end

end

@testset "lines intersects" begin
    a = Line(Point(0.0, 0.0), Point(4.0, 1.0))
    b = Line(Point(0.0, 0.25), Point(3.0, 0.25))
    c = Line(Point(0.0, 0.25), Point(0.5, 0.25))
    d = Line(Point(0.0, 0.0), Point(0.0, 4.0))
    e = Line(Point(1.0, 0.0), Point(0.0, 4.0))
    f = Line(Point(5.0, 0.0), Point(6.0, 0.0))

    @test intersects(a, b) === (true, Point(1.0, 0.25))
    @test intersects(a, c) === (false, Point(0.0, 0.0))
    @test intersects(d, d) === (false, Point(0.0, 0.0))
    found, point = intersects(d, e)
    @test found && point ≈ Point(0.0, 4.0)
    @test intersects(a, f) === (false, Point(0.0, 0.0))
end

@testset "Offsetintegers" begin
    x = 1
    @test GeometryBasicsCore.raw(x) isa Int64
    @test GeometryBasicsCore.value(x) == x

    x = ZeroIndex(1)
    @test eltype(x) == Int64

    x = OffsetInteger{0}(1)
    @test typeof(x) == OffsetInteger{0,Int64}

    x1 = OffsetInteger{0}(2)
    @test GeometryBasicsCore.pure_max(x, x1) == x1
    @test promote_rule(typeof(x), typeof(x1)) == OffsetInteger{0,Int64}
    x2 = 1
    @test promote_rule(typeof(x2), typeof(x1)) == Int64
    @test Base.to_index(x1) == 2
    @test -(x1) == OffsetInteger{0,Int64}(-2)
    @test abs(x1) == OffsetInteger{0,Int64}(2)
    @test +(x, x1) == OffsetInteger{0,Int64}(3)
    @test *(x, x1) == OffsetInteger{0,Int64}(2)
    @test -(x, x1) == OffsetInteger{0,Int64}(-1)
    #test for /
    @test div(x, x1) == OffsetInteger{0,Int64}(0)
    @test !==(x, x1)
    @test !>=(x, x1)
    @test <=(x, x1)
    @test !>(x, x1)
    @test <(x, x1)
end

@testset "view" begin
    @testset "TupleView" begin
        x = [1, 2, 3, 4, 5, 6]
        y = TupleView{2, 1}(x)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

        y = TupleView{2}(x)
        @test y == [(1, 2), (3, 4), (5, 6)]

        y = TupleView{2, 3}(x)
        @test y == [(1, 2), (4, 5)]

        y = TupleView{3, 1}(x)
        @test y == [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6)]

        y = TupleView{2, 1}(x, connect = true)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]

    end

    @testset "connected views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        x = connect(numbers, Point{2})

        @test x == Point[(1, 2), (3, 4), (5, 6)]

        line = connect(x, Line, 1)
        @test line == [Line(Point(1, 2), Point(3, 4)), Line(Point(3, 4), Point(5, 6))]

        triangles = connect(x, Triangle)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
        x = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        tetrahedra = connect(x, NSimplex{4})
        @test tetrahedra == [Tetrahedron(x[1], x[2], x[3], x[4])]
    end

    @testset "face views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})
        faces = connect([1, 2, 3], TriangleFace)
        triangles = connect(points, faces)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
        x = Point{3}(1.0)
        triangles = connect([x], [TriangleFace(1, 1, 1)])
        @test triangles == [Triangle(x, x, x)]
        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        faces = connect([1, 2, 3, 4], SimplexFace{4})
        triangles = connect(points, faces)
        @test triangles == [Tetrahedron(points...)]
    end

    @testset "reinterpret" begin
        numbers = collect(reshape(1:6, 2, 3))
        points = reinterpret(Point{2, Int}, numbers)
        @test points[1] === Point(1, 2)
        @test points[2] === Point(3, 4)
        numbers[4] = 0
        @test points[2] === Point(3, 0)
    end
end

@testset "metadata" begin
    @testset "MetaT{Point}" begin
        p = Point(1.1, 2.2)
        @test p isa AbstractVector{Float64}
        pm = MetaT(Point(1.1, 2.2); a=1, b=2)
        p1 = Point(2.2, 3.6)
        p2 = [p, p1]
        @test coordinates(p2) == p2
        @test pm.meta === (a=1, b=2)
        @test pm.main === p
        @test propertynames(pm) == (:main, :a, :b)
        @test metafree(pm) == p
        @test meta(pm) == (a = 1, b = 2)
    end

    @testset "MetaT{MultiPoint}" begin
        p = collect(Point{2, Float64}(x, x+1) for x in 1:5)
        @test p isa AbstractVector
        mpm = MetaT(MultiPoint(p); a=1, b=2)
        @test coordinates(mpm.main) == Point{2, Float64}[(x, x+1) for x in 1:5]
        @test mpm.meta === (a=1, b=2)
        @test mpm.main == p
        @test propertynames(mpm) == (:main, :a, :b)
        @test metafree(mpm) == p
        @test meta(mpm) == (a = 1, b = 2)
    end

    @testset "MetaT{LineString}" begin
        linestring = MetaT(LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)]), a = 1, b = 2)
        @test linestring isa MetaT
        @test linestring.meta === (a = 1, b = 2)
        @test propertynames(linestring) == (:main, :a, :b)
        @test metafree(linestring) == LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        @test meta(linestring) == (a = 1, b = 2)
    end

    @testset "MetaT{MultiLineString}" begin
        linestring1 = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2, Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        multilinestringmeta = MetaT(MultiLineString([linestring1, linestring2]); boundingbox = Rect(1.0, 1.0, 2.0, 2.0))
        @test multilinestringmeta isa MetaT
        @test multilinestringmeta.meta === (boundingbox = Rect(1.0, 1.0, 2.0, 2.0),)
        @test multilinestringmeta.main == multilinestring
        @test propertynames(multilinestringmeta) == (:main, :boundingbox)
        @test metafree(multilinestringmeta) == multilinestring
        @test meta(multilinestringmeta) == (boundingbox = HyperRectangle{2,Float64}([1.0, 1.0], [2.0, 2.0]),)
    end
end

@testset "tests from GeometryTypes" begin
    include("geometrytypes.jl")
end