using Test

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
