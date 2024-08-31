@testset "Meshing a single triangle sometimes returns an empty mesh" begin
    ϕ = (sqrt(5)+1)/2
    p,q,r = Point(ϕ,0,+1),Point(1,ϕ,0),Point(ϕ,0,-1)
    m = triangle_mesh(Triangle(p,q,r))
    @test m isa Mesh
    @test faces(m) == [TriangleFace(1, 2, 3)]
end

@testset "Heterogenous faces" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/142
    f = [TriangleFace(1, 2, 3), QuadFace(1, 2, 3, 4)]
    p = Point2f[(0, 1), (1, 2), (3, 4), (4, 5)]
    m = Mesh(p, f)
    @test collect(m) == [Triangle(p[1], p[2], p[3]), GeometryBasics.Quadrilateral(p[1], p[2], p[3], p[4])]
end

@testset "Heterogenous faces" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/142
    f = [TriangleFace(1, 2, 3), QuadFace(1, 2, 3, 4)]
    p = Point2f[(0, 1), (1, 2), (3, 4), (4, 5)]
    m = Mesh(p, f)
    @test collect(m) == [Triangle(p[1], p[2], p[3]), GeometryBasics.Quadrilateral(p[1], p[2], p[3], p[4])]
end

@testset "Ambiguous NgonFace constructors" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/151
    # Currently no StaticVector support
    # t = TriangleFace(SA[0.4, 0.2, 0.55])
end

@testset "Merge empty vector of meshes" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/136
    merge(Mesh[]) == Mesh(Point3f[], GLTriangleFace[])
end

@testset "Vertex Index Remapping" begin
    # Sanity Check
    # TODO: extend
    m = Mesh(
        [GeometryBasics.NormalFace(QuadFace(1, 2, 3, 4), QuadFace(1,1,1,1))],
        position = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)],
        normal = [Vec3f(0,0,1)]
    )
    
    m2 = GeometryBasics.merge_vertex_indices(m)

    @test faces(m) isa AbstractVector{<: GeometryBasics.MultiFace}
    @test propertynames(eltype(faces(m))) == keys(GeometryBasics.vertex_attributes(m))
    @test isempty(m.views)

    @test faces(m2) isa AbstractVector{<: QuadFace}
    @test coordinates(m2) == coordinates(m)
    @test normals(m2) != normals(m)
    @test normals(m2) == [only(normals(m)) for _ in 1:4]
    @test isempty(m2.views)
end