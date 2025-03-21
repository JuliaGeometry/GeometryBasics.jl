using Test
using GeometryBasics: to_ndim

@testset "Construction and Conversion" begin
    for VT in [Point, Vec]
        for T in [Int32, Float32, Float64, UInt16, BigFloat]
            p = VT{3, T}(1,2,3)
            @test p[1] == T(1)
            @test p[2] == T(2)
            @test p[3] == T(3)
        end

        for VT2 in [Point, Vec]
            @test VT{2, Float32}(VT2{3, Float32}(1,2,3)) == VT{2, Float32}(1,2)
            @test VT{2, Float32}(VT2{3, Float64}(1,2,3)) == VT{2, Float32}(1,2)
        end
        @test VT{2, Float32}(Float32[1,2,3]) == VT{2, Float32}(1,2)
        @test VT{2, Float32}([1,2,3]) == VT{2, Float32}(1,2)
    end

    @test convert(Point, (2, 3)) === Point(2, 3)
    @test convert(Point, (2.0, 3)) === Point(2.0, 3.0)

    @test to_ndim(Point3f, Vec2i(1,2), 0) == Point3f(1,2,0)
    @test to_ndim(Vec4i, (1f0, 2), 0) == Vec4i(1,2,0,0)
    @test to_ndim(NTuple{2, Float64}, Point3f(1,2,3), 0) == (1.0, 2.0)
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
                @test T1(1,2,3) .* T2(1,2,3) == T(1, 4, 9)
                # TODO: repair this:
                # @test foo.(T1(1,2,3), [T2(1,1,1), T2(2,2,2)]) == [T(1,2,3), T(2,4,6)]
                # @test foo.([T2(1,1,1), T2(2,2,2)], T1(1,2,3)) == [T(0, -3, -8), T(3, 0, -5)]
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

    @testset "chained/nested" begin
        for T1 in (Vec, Point, tuple)
            for T2 in (Vec, Point, tuple)
                T1 == tuple && T2 == tuple && continue
                T = ifelse(T1 == Point, Point, ifelse(T2 == Point, Point, Vec))
                @test T(-6, -4, 4) == T1(1,2,3) .+ T2(1, 0, 1) .- T2(3, 2, 1) .* T1(2,2,0) .- T2(2,2,0)
                @test T(-15, -5) == foo.(T1(1,2), T1(-1, 0) .+ foo.(T1(1,1), T2(2,2)))
            end
        end
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
        # @test M3[2, Vec(1,2)] == Mat{1, 2}(M3[2,1], M3[2,2])
        @test M3[2, Vec(1,2)] == Vec2(M3[2,1], M3[2,2])

        for x in (2, Vec(1,2), Vec(1,1,2,2))
            for y in (2, Vec(1,2), Vec(1,1,2,2))
                x isa Real && y isa Real && continue
                if length(x) == 1
                    @test M3[x, y] == Vec{length(y)}((M3[i, j] for j in y for i in x)...)
                elseif length(y) == 1
                    @test M3[x, y] == Vec{length(x)}((M3[i, j] for j in y for i in x)...)
                else
                    @test M3[x, y] == Mat{length(x), length(y)}((M3[i, j] for j in y for i in x)...)
                end
                @test_throws BoundsError M3[x .- 2, y]
                @test_throws BoundsError M3[x, y .+ 2]
                @test_throws BoundsError M3[x .+ 2, y .- 2]
            end
        end

    end

    for N in 1:4
        @testset "math $N x $N" begin
            bm = rand(N, N) + I
            sm = Mat{N, N}(bm)
            bv = rand(N)
            sv = Vec{N, Float64}(bv)
            @test bm == Matrix(sm)
            @test det(bm) ≈ det(sm)
            @test inv(bm) ≈ Matrix(inv(sm))
            @test collect(transpose(bm)) ≈ Matrix(transpose(sm))
            @test bm * bv ≈ collect(sm * sv)
        end
    end
end