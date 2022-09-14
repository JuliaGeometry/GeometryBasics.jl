using Test, GeometryBasicsCore

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
