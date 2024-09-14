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
        position = GeometryBasics.FaceView(Point2f[(0, 0), (1, 0), (1, 1), (0, 1)], [QuadFace(1,2,3,4)]),
        normal   = GeometryBasics.FaceView([Vec3f(0,0,1)], [QuadFace(1)])
    )
    
    m2 = GeometryBasics.clear_faceviews(m)

    @test faces(m) == [QuadFace(1,2,3,4)]
    @test coordinates(m) == Point2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    @test normals(m) == GeometryBasics.FaceView([Vec3f(0,0,1)], [QuadFace(1)])
    @test isempty(m.views)

    @test faces(m2)  == [QuadFace(1,2,3,4)]
    @test coordinates(m2) == coordinates(m)
    @test normals(m2) != normals(m)
    @test normals(m2) == [only(values(normals(m))) for _ in 1:4]
    @test isempty(m2.views)
end

@testset "complex merge" begin
    rects = [Rect3f(Point3f(x, y, z), Vec3f(0.5)) for x in -1:1 for y in -1:1 for z in -1:1]
    direct_meshes = map(rects) do r
        GeometryBasics.Mesh(coordinates(r), faces(r), normal = normals(r))
    end
    dm = merge(direct_meshes)
    
    @test GeometryBasics.facetype(dm) == QuadFace{Int64}
    @test length(faces(dm)) == 27 * 6 # 27 rects, 6 quad faces
    @test length(normals(dm)) == 27 * 6
    @test length(coordinates(dm)) == 27 * 8
    @test normals(dm) isa GeometryBasics.FaceView
    @test coordinates(dm) isa Vector
    @test !allunique([idx for f in faces(dm) for idx in f])
    @test !allunique([idx for f in faces(dm.normal) for idx in f])
    
    indirect_meshes = map(rects) do r
        GeometryBasics.mesh(coordinates(r), faces(r), normal = normals(r), facetype = QuadFace{Int64})
    end
    im = merge(indirect_meshes)

    @test im == dm

    converted_meshes = map(rects) do r
        m = GeometryBasics.Mesh(coordinates(r), faces(r), normal = normals(r))
        GeometryBasics.clear_faceviews(m)
    end
    cm = merge(converted_meshes)

    @test GeometryBasics.facetype(cm) == QuadFace{Int64}
    @test length(faces(cm)) == 27 * 6 # 27 rects, 6 quad faces
    @test length(normals(cm)) == 27 * 6 * 4 # duplicate 4x across face
    @test length(coordinates(cm)) == 27 * 8 * 3 # duplicate 3x across shared vertex
    @test normals(cm) isa Vector
    @test coordinates(cm) isa Vector
    @test allunique([idx for f in faces(cm) for idx in f])
    

    mixed_meshes = map(direct_meshes, converted_meshes) do dm, cm
        rand() > 0.5 ? dm : cm 
    end
    mm = merge(mixed_meshes)

    @test mm == cm
end

@testset "Mesh Constructor" begin
    ps = rand(Point2f, 10)
    ns = rand(Vec3f, 10)
    fs = GLTriangleFace[(1,2,3), (3,4,5), (5,6,7), (8,9,10)]

    @testset "Extracting faces from position FaceView" begin
        # can't extract from array
        @test_throws TypeError Mesh(position = ps, normal = ns)
    
        m = Mesh(position = FaceView(ps, fs), normal = ns)
        @test coordinates(m) == ps
        @test normals(m) == ns
        @test faces(m) == fs
    end

    @testset "Verifaction" begin
        # enough vertices present
        @test_throws ErrorException Mesh(rand(Point2f, 7), fs)
        m = Mesh(rand(Point2f, 12), fs)
        @test length(m.position) == 12
        @test length(m.faces) == 4
        
        @test_throws ErrorException Mesh(ps, fs, normal = rand(Vec3f, 8))
        m = Mesh(ps, fs, normal = rand(Vec3f, 12))
        @test length(m.position) == 10
        @test length(m.normal) == 12
        @test length(m.faces) == 4

        # valid FaceView (enough faces, vertices, matching dims)
        @test_throws ErrorException Mesh(ps, fs, normal = FaceView(ns, GLTriangleFace[]))
        @test_throws ErrorException Mesh(ps, fs, normal = FaceView(Vec3f[], fs))
        @test_throws ErrorException Mesh(ps, fs, normal = FaceView(ns, QuadFace{Int}.(1:4)))
        m = Mesh(ps, fs, normal = FaceView(rand(Vec3f, 9), TriangleFace{Int64}.(1:2:8)))
        @test length(m.position) == 10
        @test length(values(m.normal)) == 9
        @test length(faces(m.normal)) == 4
        @test length(m.faces) == 4

        msg = "`normals` as a vertex attribute name has been deprecated in favor of `normal` to bring it in line with mesh.position and mesh.uv"
        @test_warn msg Mesh(ps, fs, normals = ns)
    end
end

@testset "mesh() constructors" begin
    r = Rect3d(Point3d(-1), Vec3d(2))

    @testset "prerequisites" begin
        ps = coordinates(r)
        @test length(ps) == 8
        @test ps isa Vector{Point3d}
        ns = normals(r)
        @test length(ns) == 6
        @test ns isa GeometryBasics.FaceView{Vec3f, Vector{Vec3f}, Vector{QuadFace{Int64}}}
        uvs = texturecoordinates(r)
        @test length(uvs) == 8
        @test_broken uvs isa Vector{Vec2f}
        fs = faces(r)
        @test length(fs) == 6
        @test fs isa Vector{QuadFace{Int64}}
    end

    @testset "normal_mesh()" begin
        m = normal_mesh(r, pointtype = Point3f, normaltype = Vec3f)
        m = GeometryBasics.clear_faceviews(m)
        
        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 24
        @test GeometryBasics.pointtype(m) == Point3f

        @test hasproperty(m, :normal)
        @test normals(m) isa Vector{Vec3f}
        @test length(normals(m)) == 24

        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing

        @test faces(m) isa Vector{GLTriangleFace}
        @test length(faces(m)) == 12
        @test GeometryBasics.facetype(m) == GLTriangleFace
    end

    @testset "normal_uv_mesh()" begin
        m = uv_normal_mesh(
            r, pointtype = Point3d, normaltype = Vec3d, 
            uvtype = Vec3d, facetype = QuadFace{Int32}
        )

        @test hasproperty(m, :position)
        @test coordinates(m) isa Vector{Point3d}
        @test length(coordinates(m)) == 8
        @test GeometryBasics.pointtype(m) == Point3d

        @test hasproperty(m, :normal)
        @test normals(m) isa GeometryBasics.FaceView{Vec3d, Vector{Vec3d}, Vector{QuadFace{Int32}}}
        @test length(normals(m)) == 6

        @test hasproperty(m, :uv)
        @test texturecoordinates(m) isa Vector{Vec3d}
        @test length(texturecoordinates(m)) == 8

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

        @test !hasproperty(m, :normal)
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

        @test !hasproperty(m, :normal)
        @test normals(m) === nothing

        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing

        @test faces(m) isa Vector{GLTriangleFace}
        @test length(faces(m)) == 12
        @test GeometryBasics.facetype(m) == GLTriangleFace
    end

    @testset "mesh(mesh)" begin
        m = GeometryBasics.mesh(r, pointtype = Point3f, normal = normals(r), facetype = QuadFace{Int64})

        # Should be hit by normal_mesh as well...
        @test coordinates(m) isa Vector{Point3f}
        @test length(coordinates(m)) == 8
        @test normals(m) isa GeometryBasics.FaceView{Vec3f, Vector{Vec3f}, Vector{QuadFace{Int64}}}
        @test length(normals(m)) == 6
        @test !hasproperty(m, :uv)
        @test texturecoordinates(m) === nothing
        @test faces(m) isa Vector{QuadFace{Int64}}
        @test length(faces(m)) == 6

        # Shoudl throw because uv's don't match length(position) or have faces
        @test_throws ErrorException GeometryBasics.mesh(m, uv = Vec3f[])

        # remap vertex attributes to merge faceviews into one face array
        m2 = GeometryBasics.clear_faceviews(m)

        @test coordinates(m2) isa Vector{Point3f}
        @test length(coordinates(m2)) == 24
        @test normals(m2) isa Vector{Vec3f}
        @test length(normals(m2)) == 24
        @test !hasproperty(m2, :uv)
        @test texturecoordinates(m2) === nothing
        @test faces(m2) isa Vector{QuadFace{Int64}}
        @test length(faces(m2)) == 6

        # convert face type and add uvs
        m2 = GeometryBasics.mesh(m2, facetype = GLTriangleFace, uv = decompose(Point3f, m2))

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