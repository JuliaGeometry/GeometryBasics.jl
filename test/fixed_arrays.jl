using Test

@testset "conversion" begin
    @test convert(Point, (2, 3)) === Point(2, 3)
    @test convert(Point, (2.0, 3)) === Point(2.0, 3.0)
end

@testset "broadcast" begin
    x = [Point(2, 3), Point(7, 3)]

    @test [Point(4, 9), Point(14, 9)] == x .* (Point(2, 3),)
    @test [Point(4, 6), Point(9, 6)] == x .+ (Point(2, 3),)
    @test [Point(0, 0), Point(5, 0)] == x .- (Point(2, 3),)
end
