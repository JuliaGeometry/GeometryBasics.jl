using Test

@testset "conversion" begin
    @test convert(Point, (2, 3)) === Point(2, 3)
    @test convert(Point, (2.0, 3)) === Point(2.0, 3.0)
end

@testset "broadcast" begin
    foo(a, b) = (a + b) * (a - b)
    M1 = Mat{2, 2}(1,2,3,4)
    M2 = Mat{2, 2}(2,2,1,1)

    @testset "with similar" begin
        for T1 in (Vec, Point, tuple)
            for T2 in (Vec, Point, tuple)
                T1 == tuple && T2 == tuple && continue
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test T(2, 2, 4) == T1(1,2,3) .+ T2(1, 0, 1)
                @test T(foo.((1,2,3), (1, 0, 1))) == foo.(T1(1,2,3), T2(1, 0, 1))
            end
        end

        @test Mat{2, 2}(3,4,4,5) == M1 .+ M2
        @test Mat{2, 2}(foo.(values(M1), values(M2))) == foo.(M1, M2)
    end

    @testset "with const" begin
        for T in (Vec, Point)
            @test T(-4, -3) == T(1,2) .- 5
            @test T(foo.((1,2), 5)) == foo.(T(1,2), 5)
        end

        @test Mat{2, 2}(2,3,4,5) == M1 .+ 1
        @test Mat{2, 2}(foo.(values(M1), 1)) == foo.(M1, 1)
    end

    # TODO: Should this work with tuple somehow?
    @testset "with Vector (acting as const)" begin
        for T1 in (Vec, Point)
            x = [T1(2, 3), T1(7, 3)]
            for T2 in (Vec, Point)
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test [T(4, 9), T(14, 9)] == x .* T2(2, 3)
                @test [T(foo(v, T2(3, -1))) for v in x] == foo.(x, T2(3, -1))
            end
        end

        x = [M1, M2]
        @test [M1 + M2, M2 + M2] == x .+ M2
        @test [foo(M1, M2), foo(M2, M2)] == foo.(x, M2)

        # maybe bad...?
        @test [Vec(2, 3), Vec(3, 4)] == [1, 2] .+ Vec(1, 2)
    end

    @testset "chained/nested" begin
        for T1 in (Vec, Point, tuple)
            for T2 in (Vec, Point, tuple)
                T1 == tuple && T2 == tuple && continue
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test T(-6, -4, 4) == T1(1,2,3) .+ T2(1, 0, 1) .- T2(3, 2, 1) .* T1(2,2,0) .- T2(2,2,0)
                @test T(-15, -5) == foo.(T1(1,2), T1(-1, 0) .+ foo.(T1(1,1), T2(2,2)))
            end
        end

        @test [Point(1, 3), Point(3, 1)] == [Vec(1,2), Vec(2, 1)] .+ [Vec(1,2), Vec(2, 1)] .+ Point(-1, -1)
        @test [Vec(1, 3), Point(3, 1)] == [Vec(1,2), Vec(2, 1)] .+ Vec(-1, -1) .+ [Vec(1,2), Point(2, 1)]
        @test [Vec(-1, -3), Point(-3, -1)] == foo.([Vec(1,2), Vec(2, 1)] .+ Vec(-1, -1), [Vec(1,2), Point(2, 1)])

        x = [M1, M2]
        @test [M1 * M1 + M2, M1 * M2 + M2] == M1 .* x .+ M2
        @test [M1 + M2 * M1, M2 + M2 * M2] == x .+ M2 .* x
        @test [foo(M1+M2, M1), foo(M2+M2, M2)] == foo.(x .+ M2, x)
    end

    @testset "Longer functions" begin
        foo2(a, b, c, d) = (a + b) * (c + d)
        for T1 in (Vec, Point, tuple)
            for T2 in (Vec, Point, tuple)
                T1 == tuple && T2 == tuple && continue
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test T(foo2.((1,2,3), (1, 0, 1), (3, 2, 1), (2,2,0))) == 
                    foo2.(T1(1,2,3), T2(1, 0, 1), T2(3, 2, 1), T2(2,2,0))
            end
        end

        @test Mat{2, 2}(foo2.(values.((M1, M1, M2, M2))...)) == foo2.(M1, M1, M2, M2)
    end
end

@testset "finite, nan, inf tests" begin
    for T in (Vec, Point)
        @testset "$T" begin
            nan_point = T(Float64.((1.0, 2.0, 3.0, NaN)))
            inf_point = T(Float64.((1.0, 2.0, Inf, 4.0)))
            @test isinf(inf_point)
            @test !isinf(nan_point)
            @test isnan(nan_point)
            @test !isnan(inf_point)
            @test !isfinite(nan_point)
            @test !isfinite(inf_point)
            @test !isfinite(nan_point + inf_point)
            @test isfinite(T(1, 2, 3))
        end
    end
end
