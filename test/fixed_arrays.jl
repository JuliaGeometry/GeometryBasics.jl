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
            x = [T1(2, 3, 1), T1(7, 3, 2)]
            for T2 in (Vec, Point)
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test [T(4, 9, 4), T(14, 9, 8)] == x .* T2(2, 3, 4)
                @test [T(foo(v, T2(3, -1, 1))) for v in x] == foo.(x, T2(3, -1, 1))
            end
        end

        x = [M1, M2, M1]
        @test [M1 + M2, M2 + M2, M1 + M2] == x .+ M2
        @test [foo(M1, M2), foo(M2, M2), foo(M1, M2)] == foo.(x, M2)

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

@testset "Mat" begin
    M3 = Mat3(1,2,3, 4,5,6, 7,8,9)
    @test M3 isa Mat{3,3,Int,9}
    
    @testset "indexing" begin
        for i in 1:9
            @test getindex(M3, i) == i
        end
        @test_throws BoundsError getindex(M3, 0)
        @test_throws BoundsError getindex(M3, 10)

        # Sanity check for loop
        @test M3[2, Vec(1,2)] == Mat{1, 2}(M3[2,1], M3[2,2])

        for x in (2, Vec(1,2), Vec(1,1,2,2))
            for y in (2, Vec(1,2), Vec(1,1,2,2))
                x isa Real && y isa Real && continue
                @test M3[x, y] == Mat{length(x), length(y)}((M3[i, j] for j in y for i in x)...)
                @test_throws BoundsError M3[x .- 2, y]
                @test_throws BoundsError M3[x, y .+ 2]
                @test_throws BoundsError M3[x .+ 2, y .- 2]
            end
        end

    end
end