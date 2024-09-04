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
        [GeometryBasics.MultiFace(position = QuadFace(1, 2, 3, 4), normal = QuadFace(1,1,1,1))],
        position = Point2f[(0, 0), (1, 0), (1, 1), (0, 1)],
        normals = [Vec3f(0,0,1)]
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

@testset "mesh() constructors" begin
    r = Rect3d(Point3d(-1), Vec3d(2))

    @testset "prerequisites" begin
        ps = collect(coordinates(r))
        @test length(ps) == 8
        @test ps isa Vector{Point3d}
        ns = collect(normals(r))
        @test length(ns) == 6
        @test ns isa Vector{Vec3f}
        uvs = collect(texturecoordinates(r))
        @test length(uvs) == 8
        @test_broken uvs isa Vector{Vec2f}
        fs = collect(faces(r))
        @test length(fs) == 6
        @test fs isa Vector{GeometryBasics.NormalUVFace{4, Int64, QuadFace{Int64}}}
    end

    @testset "normal_mesh()" begin
        # TODO: simplify?
        FT = GeometryBasics.NormalFace{4, Int64, QuadFace{Int64}}
        m = normal_mesh(r, pointtype = Point3f, normaltype = Vec3f, facetype = FT)
        
        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 8
        @test GeometryBasics.pointtype(m) == Point3f

        @test hasproperty(m, :normals)
        @test normals(m) isa Vector{Vec3f}
        @test length(normals(m)) == 6

        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing

        @test faces(m) isa Vector{FT}
        @test length(faces(m)) == 6
        @test GeometryBasics.facetype(m) == FT
    end

    @testset "normal_uv_mesh()" begin
        m = uv_normal_mesh(
            r, pointtype = Point3d, normaltype = Vec3d, 
            uvtype = Vec3d, facetype = QuadFace{Int32}
        )

        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3d}
        @test length(coordinates(m)) == 24
        @test GeometryBasics.pointtype(m) == Point3d

        @test hasproperty(m, :normals)
        @test normals(m) isa Vector{Vec3d}
        @test length(normals(m)) == 24

        @test hasproperty(m, :uv)
        @test texturecoordinates(m) isa Vector{Vec3d}
        @test length(texturecoordinates(m)) == 24

        @test faces(m) isa Vector{QuadFace{Int32}}
        @test length(faces(m)) == 6
        @test GeometryBasics.facetype(m) == QuadFace{Int32}
    end

    @testset "uv_mesh()" begin
        m = uv_mesh(
            r, pointtype = Point3f, uvtype = Vec3f, facetype = GLTriangleFace
        )

        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 8
        @test GeometryBasics.pointtype(m) == Point3f

        @test !hasproperty(m, :normals)
        @test normals(m) === nothing

        @test hasproperty(m, :uv)
        @test texturecoordinates(m) isa Vector{Vec3f}
        @test length(texturecoordinates(m)) == 8

        @test faces(m) isa Vector{GLTriangleFace}
        @test length(faces(m)) == 12
        @test GeometryBasics.facetype(m) == GLTriangleFace
    end

    @testset "triangle_mesh()" begin
        m = triangle_mesh(r, pointtype = Point3f)

        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 8
        @test GeometryBasics.pointtype(m) == Point3f

        @test !hasproperty(m, :normals)
        @test normals(m) === nothing

        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing

        @test faces(m) isa Vector{GLTriangleFace}
        @test length(faces(m)) == 12
        @test GeometryBasics.facetype(m) == GLTriangleFace
    end

    @testset "mesh(mesh)" begin
        FT = GeometryBasics.NormalFace{4, Int64, QuadFace{Int64}}
        m = GeometryBasics.mesh(r, pointtype = Point3f, normals = normals(r), facetype = FT)

        # Should be hit by normal_mesh as well...
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 8
        @test normals(m) isa Vector{Vec3f}
        @test length(normals(m)) == 6
        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing
        @test faces(m) isa Vector{FT}
        @test length(faces(m)) == 6

        # shouldn't be able to add vertex attributes when they aren't synchronized
        @test_throws ErrorException GeometryBasics.mesh(m, uv = Vec3f[])

        # Can convert face type - should remap
        m2 = GeometryBasics.mesh(m, facetype = QuadFace{Int32})

        @test coordinates(m2) isa Vector{Point3f}
        @test length(coordinates(m2)) == 24
        @test normals(m2) isa Vector{Vec3f}
        @test length(normals(m2)) == 24
        @test !hasproperty(m2, :uv)
        @test texturecoordinates(m2) === nothing
        @test faces(m2) isa Vector{QuadFace{Int32}}
        @test length(faces(m2)) == 6

        # remap + decompose face and then adding vertex attribute is fine
        m2 = GeometryBasics.mesh(m, facetype = GLTriangleFace)
        # should be same length as other vertices...
        @test_throws ErrorException GeometryBasics.mesh(m2, uv = Vec3f[])
        m2 = GeometryBasics.mesh(m2, uv = decompose(Point3f, m2))

        @test coordinates(m2) isa Vector{Point3f}
        @test length(coordinates(m2)) == 24
        @test normals(m2) isa Vector{Vec3f}
        @test length(normals(m2)) == 24
        @test texturecoordinates(m2) isa Vector{Point3f}
        @test length(texturecoordinates(m2)) == 24
        @test faces(m2) isa Vector{GLTriangleFace}
        @test length(faces(m2)) == 12


    end
end