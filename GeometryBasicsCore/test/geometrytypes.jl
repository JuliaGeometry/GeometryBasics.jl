using Test, GeometryBasicsCore

@testset "Cylinder" begin
    @testset "constructors" begin
        o, extr, r = Point2f(1, 2), Point2f(3, 4), 5.0f0
        s = Cylinder(o, extr, r)
        @test typeof(s) == Cylinder{2,Float32}
        @test typeof(s) == Cylinder2{Float32}
        @test origin(s) == o
        @test extremity(s) == extr
        @test radius(s) == r
        #@test abs(height(s)- norm([1,2]-[3,4]))<1e-5
        h = norm(o - extr)
        @test isapprox(height(s), h)
        #@test norm(direction(s) - Point{2,Float32}([2,2]./norm([1,2]-[3,4])))<1e-5
        @test isapprox(direction(s), Point2f(2, 2) ./ h)
        v1 = rand(Point{3,Float64})
        v2 = rand(Point{3,Float64})
        R = rand()
        s = Cylinder(v1, v2, R)
        @test typeof(s) == Cylinder{3,Float64}
        @test typeof(s) == Cylinder3{Float64}
        @test origin(s) == v1
        @test extremity(s) == v2
        @test radius(s) == R
        @test height(s) == norm(v2 - v1)
        #@test norm(direction(s) - Point{3,Float64}((v2-v1)./norm(v2-v1)))<1e-10
        @test isapprox(direction(s), (v2 - v1) ./ norm(v2 .- v1))
    end
end

@testset "Rectangles" begin
    rect = Rect2f(0, 7, 20, 3)
    @test (rect + 4) == Rect2f(4, 11, 20, 3)
    @test (rect + Vec(2, -2)) == Rect2f(2, 5, 20, 3)

    @test (rect - 4) == Rect2f(-4, 3, 20, 3)
    @test (rect - Vec(2, -2)) == Rect2f(-2, 9, 20, 3)

    base = Vec3f(1, 2, 3)
    wxyz = Vec3f(-2, 4, 2)
    rect = Rect3f(base, wxyz)
    @test (rect + 4) == Rect3f(base .+ 4, wxyz)
    @test (rect + Vec(2, -2, 3)) == Rect3f(base .+ Vec(2, -2, 3), wxyz)

    @test (rect - 4) == Rect3f(base .- 4, wxyz)
    @test (rect - Vec(2, -2, 7)) == Rect3f(base .- Vec(2, -2, 7), wxyz)

    rect = Rect2f(0, 7, 20, 3)
    @test (rect * 4) == Rect2f(0, 7 * 4, 20 * 4, 3 * 4)
    @test (rect * Vec(2, -2)) == Rect2f(0, -7 * 2, 20 * 2, -3 * 2)

    base = Vec3f(1, 2, 3)
    wxyz = Vec3f(-2, 4, 2)
    rect = Rect3f(base, wxyz)
    @test (rect * 4) == Rect3f(base .* 4, wxyz .* 4)
    @test (rect * Vec(2, -2, 3)) == Rect3f(base .* Vec(2, -2, 3), wxyz .* Vec(2, -2, 3))

    rect1 = Rect(Vec(0.0, 0.0), Vec(1.0, 2.0))
    rect2 = Rect(0.0, 0.0, 1.0, 2.0)
    @test rect1 isa HyperRectangle{2,Float64}
    @test rect1 == rect2

    split1, split2 = GeometryBasicsCore.split(rect1, 2, 1)
    @test widths(split1) == widths(split2)
    @test origin(split1) == Vec(0, 0)
    @test origin(split2) == Vec(0, 1)
    @test in(split1, rect1)
    @test !in(rect1, split1)

    prim = Rect(0.0, 0.0, 1.0, 1.0)
    @test length(prim) == 2

    @test width(prim) == 1.0
    @test height(prim) == 1.0

    b1 = Rect2(0.0, 0.0, 2.0, 2.0)
    b2 = Rect2(0, 0, 2, 2)
    @test isequal(b1, b2)

    pt = Point(1.0, 1.0)
    b1 = Rect(0.0, 0.0, 1.0, 1.0)
    @test in(pt, b1)

    rect = Rect(0.0, 0.0, 1.0, 1.0)
    @test GeometryBasicsCore.positive_widths(rect) isa HyperRectangle{2,Float64}

    h1 = Rect(0.0, 0.0, 1.0, 1.0)
    h2 = Rect(1.0, 1.0, 2.0, 2.0)
    @test union(h1, h2) isa HyperRectangle{2,Float64}
    @test GeometryBasicsCore.diff(h1, h2) == h1
    @test GeometryBasicsCore.intersect(h1, h2) isa HyperRectangle{2,Float64}

    b = Rect(0.0, 0.0, 1.0, 1.0)
    v = Vec(1, 2)
    @test update(b, v) isa HyperRectangle{2,Float64}
    v = Vec(1.0, 2.0)
    @test update(b, v) isa HyperRectangle{2,Float64}

    p = Vec(5.0, 4.0)
    rect = Rect(0.0, 0.0, 1.0, 1.0)
    @test min_dist_dim(rect, p, 1) == 4.0
    @test min_dist_dim(rect, p, 2) == 3.0
    @test max_dist_dim(rect, p, 1) == 5.0
    @test max_dist_dim(rect, p, 2) == 4.0

    rect1 = Rect(0.0, 0.0, 1.0, 1.0)
    rect2 = Rect(3.0, 1.0, 4.0, 2.0)
    @test min_dist_dim(rect1, rect2, 1) == 2.0
    @test min_dist_dim(rect1, rect2, 2) == 0.0
    @test max_dist_dim(rect1, rect2, 1) == 7.0
    @test max_dist_dim(rect1, rect2, 2) == 3.0

    @test !before(rect1, rect2)
    rect1 = Rect(0.0, 0.0, 1.0, 1.0)
    rect2 = Rect(3.0, 2.0, 4.0, 2.0)
    @test before(rect1, rect2)

    @test !meets(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 4.0, 2.0)
    @test meets(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !overlaps(rect1, rect2)
    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(1.5, 1.5, 2.0, 2.0)
    @test overlaps(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !GeometryBasicsCore.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 1.5, 1.5)
    @test !GeometryBasicsCore.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 3.0, 3.0)
    @test GeometryBasicsCore.starts(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 4.0, 4.0)
    @test during(rect1, rect2)
    rect1 = Rect(0.0, 0.0, 2.0, 3.0)
    rect2 = Rect(1.0, 1.0, 4.0, 2.0)
    @test !during(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 4.0, 4.0)
    @test !finishes(rect1, rect2)
    rect1 = Rect(1.0, 0.0, 1.0, 1.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !finishes(rect1, rect2)
    rect1 = Rect(1.0, 1.0, 1.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 3.0)
    @test finishes(rect1, rect2)

end
