@testset "Arithmetics" begin
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
end

@testset "Matrix Multiplications" begin
    r = Rect2f(-1, -2, 4, 3)

    # TODO: this seems quite dangerous: We pad points with ones which makes
    #       sense for translations if we go to D+1, but is nonsense if we
    #       go higher dimensions than that.
    M = rand(Mat4f)
    ps = Point2f[M * Point(p..., 1, 1) for p in coordinates(r)]
    @test Rect2f(ps) ≈ M * r

    M = Mat2f(0.5, -0.3, 0.7, 1.5)
    ps = Point2f[M * p for p in coordinates(r)]
    @test Rect2f(ps) ≈ M * r

    r = Rect3f(-1, -2, -3, 2, 4, 1)
    M = rand(Mat4f)
    ps = Point3f[M * Point(p..., 1) for p in coordinates(r)]
    @test Rect3f(ps) ≈ M * r
end

@testset "Constructors & Getters" begin
    rect1 = Rect(Vec(0.0, 0.0), Vec(1.0, 2.0))
    rect2 = Rect(0.0, 0.0, 1.0, 2.0)
    @test rect1 isa GeometryBasics.HyperRectangle{2,Float64}
    @test rect1 == rect2

    prim = Rect(0.0, 0.0, 1.0, 1.0)
    @test length(prim) == 2

    @test width(prim) == 1.0
    @test height(prim) == 1.0

    rect = Rect(0.0, 0.0, 1.0, 1.0)
    @test GeometryBasics.positive_widths(rect) == rect
    @test GeometryBasics.positive_widths(Rect2f(0, 0, -1, -2)) == Rect2f(-1, -2, 1, 2)

    rect1 = @inferred Rect(1, 2, 3, 4, 5, 6, 7, 8)
    rect2 = Rect(Vec(1, 2, 3, 4), Vec(5, 6, 7, 8))
    @test rect1 == rect2
end

@testset "checks" begin
    @testset "isempty" begin
        @test isempty(Rect2f())
        @test isempty(Rect3f())
        @test isempty(Rect2i())
        @test isempty(Rect(0, 0, 0, 0))
        @test isempty(Rect(0, 0, 0, 0, 0, 0))
        @test isempty(Rect(0, 0, 0, 0, 0, 0))
    end

    @testset "equality" begin
        b1 = Rect2(0.0, 0.0, 2.0, 2.0)
        b2 = Rect2(0, 0, 2, 2)
        @test isequal(b1, b2)
        @test !isequal(b1, Rect2f(0.1, 0, 2, 2))
    end

    @testset "in" begin
        pt = Point(1.0, 1.0)
        b1 = Rect(0.0, 0.0, 1.0, 1.0)
        @test in(pt, b1)
        @test in(Rect2f(0, 0, 1, 1), Rect2f(-1, -1, 3, 3))
        @test in(Rect2f(0, 0, 1, 1), Rect2f(0, 0, 3, 3))
        @test !in(Rect2f(0, 0, 1, 1), Rect2f(0.1, 0, 3, 3))
    end

    @testset "overlaps" begin
        @test !overlaps(Rect2f(-1, -1, 2, 2), Rect2f(-2, -2, 1, 1))
        @test !overlaps(Rect2f(-1, -1, 2, 2), Rect2f(-2, -2, 1, 1.2))
        @test overlaps(Rect2f(-1, -1, 2, 2), Rect2f(-2, -2, 1.2, 1.2))
        @test overlaps(Rect2f(-1, -1, 2, 2), Rect2f(-2, -2, 4, 4))
        @test overlaps(Rect2f(-1, -1, 2, 2), Rect2f(-0.5, -0.5, 1, 1))
        @test overlaps(Rect2f(-1, -1, 2, 2), Rect2f(0.5, 0.5, 1, 1))
        @test !overlaps(Rect2f(-1, -1, 2, 2), Rect2f(1, 1, 1, 1))
        @test !overlaps(Rect2f(-1, -1, 2, 2), Rect2f(2, 2, 1, 1))

        @test !overlaps(Rect2f(-2, -0.1, 0.5, 0.2), Rect2f(-1, -1, 2, 2))
        @test !overlaps(Rect2f(-2, -0.1, 1, 0.2), Rect2f(-1, -1, 2, 2))
        @test overlaps(Rect2f(-2, -0.1, 1.1, 0.2), Rect2f(-1, -1, 2, 2))
        @test overlaps(Rect2f(-2, -0.1, 3, 0.2), Rect2f(-1, -1, 2, 2))
        @test overlaps(Rect2f(-1, -0.1, 3, 0.2), Rect2f(-1, -1, 2, 2))
        @test overlaps(Rect2f(0, -0.1, 3, 0.2), Rect2f(-1, -1, 2, 2))
        @test !overlaps(Rect2f(1, -0.1, 3, 0.2), Rect2f(-1, -1, 2, 2))
        @test !overlaps(Rect2f(2, -0.1, 3, 0.2), Rect2f(-1, -1, 2, 2))
    end
end

@testset "Set Operations" begin
    @testset "split" begin
        rect1 = Rect(Vec(0.0, 0.0), Vec(1.0, 2.0))
        split1, split2 = GeometryBasics.split(rect1, 2, 1)
        @test widths(split1) == widths(split2)
        @test origin(split1) == Vec(0, 0)
        @test origin(split2) == Vec(0, 1)
        @test in(split1, rect1) && in(split2, rect1)
        @test !(in(rect1, split1) || in(rect1, split2))

        rect1 = Rect(Vec(0.0, 0.0, -1.0), Vec(1.0, 2.0, 1.0))
        split1, split2 = GeometryBasics.split(rect1, 1, 0.75)
        @test widths(split1) == Vec(0.75, 2, 1)
        @test widths(split2) == Vec(0.25, 2, 1)
        @test origin(split1) == Vec(0, 0, -1)
        @test origin(split2) == Vec(0.75, 0, -1)
        @test in(split1, rect1) && in(split2, rect1)
        @test !(in(rect1, split1) || in(rect1, split2))
    end

    @testset "union" begin
        h1 = Rect(0.0, 0.0, 1.0, 1.0)
        h2 = Rect(1.0, 1.0, 2.0, 2.0)
        @test union(h1, h2) == Rect2d(0, 0, 3, 3)
        @test union(Rect2f(), Rect2f(0, 0, 1, 1)) == Rect2f(0, 0, 1, 1)
        @test union(Rect2f(0, 0, 1, 1), Rect2f()) == Rect2f(0, 0, 1, 1)
        @test union(Rect3f(0,0,0,1,1,1), Rect3f(-2,-2,-2,1,1,1)) == Rect3f(-2,-2,-2,3,3,3)
    end

    @testset "intersection" begin
        h1 = Rect(0.0, 0.0, 1.0, 1.0)
        h2 = Rect(1.0, 1.0, 2.0, 2.0)
        @test intersect(h1, h2) == Rect()
        @test intersect(Rect2f(0,0,2,1), Rect2f(1,0,1,2)) == Rect2f(1, 0, 1, 1)
        @test intersect(Rect2f(0,0,2,2), Rect2f(1,1,2,2)) == Rect2f(1, 1, 1, 1)
        @test intersect(Rect2f(0,0,2,3), Rect2f(1,1,2,1)) == Rect2f(1, 1, 1, 1)
        @test intersect(Rect2f(0,0,3,3), Rect2f(1,1,1,1)) == Rect2f(1, 1, 1, 1)
    end

    @testset "diff" begin
        a = Rectf(-1, -1, -1, 2, 2, 2)
        # no intersection -> return a
        @test bbox_diff(a, Rectf(2, -1, -1, 2, 2, 2)) == a
        @test bbox_diff(a, Rectf(2, 2, 2, 3, 3, 3)) == a
        @test bbox_diff(a, Rectf(-3, -3, -3, 1, 1, 1)) == a
        @test bbox_diff(a, Rectf(-8, -4, -4, 4, 10, 8)) == a

        # intersection is a chunk < a in every dimension -> return a
        @test bbox_diff(a, Rectf(0, 0, 0, 2, 2, 2)) == a
        @test bbox_diff(a, Rectf(-2, -4, -1, 2, 3.8, 1)) == a
        @test bbox_diff(a, Rectf(-0.5, -0.5, -1, 1, 1, 2)) == a

        # intersection is matching the size of a in just one dimension -> return a
        @test bbox_diff(a, Rectf(-2, 0, 0, 4, 2, 2)) == a
        @test bbox_diff(a, Rectf(-0.5, -0.5, -2, 1, 1, 10)) == a

        # intersection is a slab matching the size of a in two dimension -> reduced
        @test bbox_diff(a, Rectf(-2, -2, 0, 4, 4, 2)) == Rectf(-1, -1, -1, 2, 2, 1)
        @test bbox_diff(a, Rectf(-5, -2, -2, 5.2, 4, 4)) ≈ Rectf(0.2, -1, -1, 0.8, 2, 2)
        # edge case: bisection
        @test bbox_diff(a, Rectf(-0.5, -2, -2, 1, 4, 4)) == a

        # intersection is a -> empty Rect
        @test bbox_diff(a, Rectf(-2, -2, -2, 4, 4, 4)) == Rect3f()
        @test bbox_diff(a, Rectf(-1, -1, -1, 2, 2, 2)) == Rect3f()
    end
end

@testset "update" begin
    b = Rect(0.0, 0.0, 1.0, 1.0)
    @test update(b, Vec(1, 2)) == Rect2d(0, 0, 1, 2)
    @test update(b, Vec(2.0, 2.0)) == Rect2d(0, 0, 2, 2)
    @test update(b, Vec(-1, 2)) == Rect2d(-1, 0, 2, 2)
end

@testset "euclidean distances" begin
    p = Vec(5.0, 4.0)
    rect = Rect(0.0, 0.0, 1.0, 1.0)
    @test min_dist_dim(rect, p, 1) == 4.0
    @test min_dist_dim(rect, p, 2) == 3.0
    @test max_dist_dim(rect, p, 1) == 5.0
    @test max_dist_dim(rect, p, 2) == 4.0
    @test minmax_dist_dim(rect, p, 1) == (4.0, 5.0)

    rect1 = Rect(0.0, 0.0, 1.0, 1.0)
    rect2 = Rect(3.0, 1.0, 4.0, 2.0)
    @test min_dist_dim(rect1, rect2, 1) == 2.0
    @test min_dist_dim(rect1, rect2, 2) == 0.0
    @test max_dist_dim(rect1, rect2, 1) == 7.0
    @test max_dist_dim(rect1, rect2, 2) == 3.0
    @test minmax_dist_dim(rect1, rect2, 1) == (2.0, 7.0)

    r = Rect2f(-1, -1, 2, 3)
    p = Point2f(1, 2) + Point2f(3, 4)
    @test min_euclidean(r, p) == 5f0
    @test max_euclidean(r, p) ≈ sqrt(5*5 + 7*7)

    r2 = Rect2f(0, 0, 2, 3)
    @test min_euclidean(r, r2) == 0f0
    @test max_euclidean(r, r2) == 5f0
    @test minmax_euclidean(r, r2) == (0f0, 5f0)
end

# Does this make sense to extend to Rects?
@testset "Allen's interval algebra" begin
    rect1 = Rect(0.0, 0.0, 1.0, 1.0)
    rect2 = Rect(3.0, 1.0, 4.0, 2.0)
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
    @test !GeometryBasics.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 1.5, 1.5)
    @test !GeometryBasics.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 3.0, 3.0)
    @test GeometryBasics.starts(rect1, rect2)

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

@testset "Utilities" begin
    # TODO: this is effectively 0-indexed... should it be?
    M = reshape(collect(11:100), 10, 9)[1:9, :]
    r = Rect2i(2, 4, 2, 4)
    @test M[r] == [53 63 73 83; 54 64 74 84]
end