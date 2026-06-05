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

    facedata = FaceView([:red, :blue], [TriangleFace(1), QuadFace(2)])
    m2 = GeometryBasics.mesh(m, color = facedata)
    m3 = expand_faceviews(m2)
    @test faces(m3) == GLTriangleFace[(1,2,3), (4,5,6), (4,6,7)]
    @test coordinates(m3) == Point2f[[0.0, 1.0], [1.0, 2.0], [3.0, 4.0], [0.0, 1.0], [1.0, 2.0], [3.0, 4.0], [4.0, 5.0]]
    @test m3.color == [:red, :red, :red, :blue, :blue, :blue, :blue]

    @test per_face([:red, :blue], f) == facedata
    @test per_face([:red, :blue], m) == facedata
    @test per_face([:red, :blue, :blue], m2) == FaceView([:red, :blue, :blue], GLTriangleFace.(1:3))
    @test per_face([:red, :blue, :blue], m3) == FaceView([:red, :blue, :blue], GLTriangleFace.(1:3))
end

@testset "Ambiguous NgonFace constructors" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/151
    # Currently no StaticVector support
    # t = TriangleFace(SA[0.4, 0.2, 0.55])
end

@testset "Merge empty vector of meshes" begin
    # https://github.com/JuliaGeometry/GeometryBasics.jl/issues/136
    merge(Mesh[]) == Mesh(Point3f[], GLTriangleFace[])
    merge([Mesh(Point3f[], GLTriangleFace[])]) == Mesh(Point3f[], GLTriangleFace[])
    merge([Mesh(Point3f[], GLTriangleFace[]), Mesh(Point3f[], GLTriangleFace[])]) == Mesh(Point3f[], GLTriangleFace[])
end

@testset "Vertex Index Remapping" begin
    # Sanity Check
    m = Mesh(
        position = GeometryBasics.FaceView(Point2f[(0, 0), (1, 0), (1, 1), (0, 1)], [QuadFace(1,2,3,4)]),
        normal   = GeometryBasics.FaceView([Vec3f(0,0,1)], [QuadFace(1)])
    )

    m2 = GeometryBasics.expand_faceviews(m)

    @test faces(m) == [QuadFace(1,2,3,4)]
    @test coordinates(m) == Point2f[(0, 0), (1, 0), (1, 1), (0, 1)]
    @test normals(m) == GeometryBasics.FaceView([Vec3f(0,0,1)], [QuadFace(1)])
    @test isempty(m.views)

    @test per_face([Vec3f(0,0,1)], m) == m.normal

    @test faces(m2)  == [QuadFace(1,2,3,4)]
    @test coordinates(m2) == coordinates(m)
    @test normals(m2) != normals(m)
    @test normals(m2) == [only(values(normals(m))) for _ in 1:4]
    @test isempty(m2.views)

    mm = MetaMesh(m, Dict(:test => 1, :a => "a"))
    mm2 = GeometryBasics.expand_faceviews(mm)
    @test mm2.meta == mm.meta
    @test Mesh(mm2) == m2
end

@testset "Duplicate face removal" begin
    fs = GLTriangleFace[(1,2,3), (2,3,4), (3,4,5), (1,2,3), (1,4,5)]
    fs = [fs; fs]
    new_fs = remove_duplicates(fs)
    @test all(f -> f in GLTriangleFace[(1,2,3), (2,3,4), (3,4,5), (1,4,5)], new_fs)

    fs = rand(QuadFace{Int32}, 4)
    new_fs = remove_duplicates([fs; fs])
    @test all(in(fs), new_fs)
end

@testset "complex merge + split" begin
    rects = [Rect3f(Point3f(x, y, z), Vec3f(0.5)) for x in -1:1 for y in -1:1 for z in -1:1]
    direct_meshes = map(rects) do r
        GeometryBasics.Mesh(coordinates(r), faces(r), normal = normals(r))
    end
    dm = merge(direct_meshes)

    @test GeometryBasics.facetype(dm) == QuadFace{Int}
    @test length(faces(dm)) == 27 * 6 # 27 rects, 6 quad faces
    @test length(normals(dm)) == 27 * 6
    @test length(coordinates(dm)) == 27 * 8
    @test normals(dm) isa GeometryBasics.FaceView
    @test coordinates(dm) isa Vector
    @test !allunique([idx for f in faces(dm) for idx in f])
    @test !allunique([idx for f in faces(dm.normal) for idx in f])

    split_meshes = split_mesh(dm)
    @test all(GeometryBasics.strictly_equal_face_vertices.(split_meshes, direct_meshes))

    indirect_meshes = map(rects) do r
        m = GeometryBasics.mesh(coordinates(r), faces(r), normal = normals(r), facetype = QuadFace{Int64})
        # Also testing merge of meshes with views
        push!(m.views, 1:length(faces(m)))
        m
    end
    im = merge(indirect_meshes)

    @test im == dm
    @test GeometryBasics.facetype(im) == QuadFace{Int64}

    split_meshes = split_mesh(im)
    @test all(GeometryBasics.strictly_equal_face_vertices.(split_meshes, indirect_meshes))

    converted_meshes = map(rects) do r
        m = GeometryBasics.Mesh(coordinates(r), faces(r), normal = normals(r))
        GeometryBasics.expand_faceviews(m)
    end
    cm = merge(converted_meshes)

    @test GeometryBasics.facetype(cm) == QuadFace{Int}
    @test length(faces(cm)) == 27 * 6 # 27 rects, 6 quad faces
    @test length(normals(cm)) == 27 * 6 * 4 # duplicate 4x across face
    @test length(coordinates(cm)) == 27 * 8 * 3 # duplicate 3x across shared vertex
    @test normals(cm) isa Vector
    @test coordinates(cm) isa Vector
    @test allunique([idx for f in faces(cm) for idx in f])

    split_meshes = split_mesh(cm)
    @test all(GeometryBasics.strictly_equal_face_vertices.(split_meshes, converted_meshes))

    mixed_meshes = map(direct_meshes, indirect_meshes, converted_meshes) do dm, im, cm
        rand((dm, im, cm)) # (with FaceView, with mesh.views & FaceView, w/o FaceView)
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

    @test_throws MethodError Mesh(Point2f[], GLTriangleFace(1,2,3))

    @testset "Verification" begin
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
        @test_logs (:warn, msg) Mesh(ps, fs, normals = ns)
    end
end

@testset "Interface" begin
    ps = rand(Point2f, 10)
    ns = rand(Vec3f, 10)
    uvs = FaceView(rand(Vec2f, 4), GLTriangleFace.(1:4))
    fs = GLTriangleFace[(1,2,3), (3,4,5), (5,6,7), (8,9,10)]

    m = Mesh(ps, fs, normal = ns, uv = uvs)

    @test vertex_attributes(m) == getfield(m, :vertex_attributes)
    @test coordinates(m) == vertex_attributes(m)[:position]
    @test normals(m) == vertex_attributes(m)[:normal]
    @test texturecoordinates(m) == vertex_attributes(m)[:uv]
    @test faces(m) == getfield(m, :faces)

    @test m.vertex_attributes == getfield(m, :vertex_attributes)
    @test m.position == vertex_attributes(m)[:position]
    @test m.normal == vertex_attributes(m)[:normal]
    @test m.uv == vertex_attributes(m)[:uv]
    @test m.faces == getfield(m, :faces)

    @test hasproperty(m, :vertex_attributes)
    @test hasproperty(m, :position)
    @test hasproperty(m, :normal)
    @test hasproperty(m, :uv)
    @test hasproperty(m, :faces)

    @test GeometryBasics.meta(m) == NamedTuple()
    @test Mesh(m) === m

    mm = MetaMesh(m, name = "test")

    @test Mesh(mm) == m
    @test haskey(mm, :name)
    @test get(mm, :name, nothing) == "test"
    @test get(() -> nothing, mm, :name) == "test"
    @test get(mm, :foo, nothing) === nothing
    @test get(() -> nothing, mm, :foo) === nothing
    @test mm[:name] == "test"
    @test !haskey(mm, :foo)
    @test get!(mm, :foo, "foo") == "foo"
    @test haskey(mm, :foo)
    @test get!(() -> "bar", mm, :bar) == "bar"
    @test haskey(mm, :bar)
    mm[:foo] = "bar"
    @test mm[:foo] == "bar"
    delete!(mm, :bar)
    @test !haskey(mm, :bar)
    @test keys(mm) == keys(getfield(mm, :meta))

    @test vertex_attributes(mm) == getfield(m, :vertex_attributes)
    @test coordinates(mm) == vertex_attributes(m)[:position]
    @test normals(mm) == vertex_attributes(m)[:normal]
    @test texturecoordinates(mm) == vertex_attributes(m)[:uv]
    @test faces(mm) == getfield(m, :faces)

    @test hasproperty(mm, :vertex_attributes)
    @test hasproperty(mm, :position)
    @test hasproperty(mm, :normal)
    @test hasproperty(mm, :uv)
    @test hasproperty(mm, :faces)
    @test propertynames(mm) == (:mesh, :meta, propertynames(m)...)

    @test mm.vertex_attributes == getfield(m, :vertex_attributes)
    @test mm.position == vertex_attributes(m)[:position]
    @test mm.normal == vertex_attributes(m)[:normal]
    @test mm.uv == vertex_attributes(m)[:uv]
    @test mm.faces == getfield(m, :faces)

    @test GeometryBasics.meta(mm) == mm.meta
    @test Mesh(mm) === mm.mesh
end

@testset "mesh() constructors" begin
    r = Rect3d(Point3d(-1), Vec3d(2))

    @testset "prerequisites" begin
        ps = coordinates(r)
        @test length(ps) == 8
        @test ps isa Vector{Point3d}
        ns = normals(r)
        @test length(ns) == 6
        @test ns isa GeometryBasics.FaceView{Vec3f, Vector{Vec3f}, Vector{QuadFace{Int}}}
        uvs = texturecoordinates(r)
        @test length(uvs) == 8
        @test_broken uvs isa Vector{Vec2f}
        fs = faces(r)
        @test length(fs) == 6
        @test fs isa Vector{QuadFace{Int}}
    end

    @testset "normal_mesh()" begin
        m = normal_mesh(r, pointtype = Point3f, normaltype = Vec3f)
        m = GeometryBasics.expand_faceviews(m)

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

        @test normal_mesh(coordinates(m), faces(m)) == m
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
        m2 = GeometryBasics.expand_faceviews(m)

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

@testset "map_coordinates" begin
    m = GeometryBasics.mesh(Rect3f(0,0,0,1,1,1))
    m2 = GeometryBasics.map_coordinates(p -> 2 * p, m)
    @test m !== m2
    @test 2 * coordinates(m) == coordinates(m2)

    m3 = GeometryBasics.map_coordinates!(p -> 0.5 * p, m2)
    @test m3 === m2
    @test coordinates(m) == coordinates(m3)

    m = MetaMesh(GeometryBasics.mesh(Rect3f(0,0,0,1,1,1)), meta = "test")
    m2 = GeometryBasics.map_coordinates(p -> 2 * p, m)
    @test m !== m2
    @test 2 * coordinates(m) == coordinates(m2)
    @test GeometryBasics.meta(m) == GeometryBasics.meta(m2)
end