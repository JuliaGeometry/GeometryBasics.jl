using Test, GeometryBasics

@testset "Cylinder" begin
    @testset "constructors" begin
        v1 = rand(Point{3,Float64})
        v2 = rand(Point{3,Float64})
        R = rand()
        s = Cylinder(v1, v2, R)
        @test typeof(s) == Cylinder{Float64}
        @test origin(s) == v1
        @test extremity(s) == v2
        @test radius(s) == R
        @test height(s) == norm(v2 - v1)
        #@test norm(direction(s) - Point{3,Float64}((v2-v1)./norm(v2-v1)))<1e-10
        @test isapprox(direction(s), (v2 - v1) ./ norm(v2 .- v1))
    end

    @testset "decompose" begin
        v1 = Point{3,Float64}(1, 2, 3)
        v2 = Point{3,Float64}(4, 5, 6)
        R = 5.0
        s = Cylinder(v1, v2, R)
        positions = Point{3,Float64}[
            (4.535533905932738, -1.5355339059327373, 3.0),
            (3.0412414523193148, 4.041241452319315, -1.0824829046386295),
            (-2.535533905932737, 5.535533905932738, 2.9999999999999996),
            (-1.0412414523193152, -0.04124145231931431, 7.0824829046386295),
            (7.535533905932738, 1.4644660940672627, 6.0),
            (6.041241452319315, 7.041241452319315, 1.9175170953613705),
            (0.46446609406726314, 8.535533905932738, 6.0),
            (1.9587585476806848, 2.9587585476806857, 10.08248290463863),
            (1, 2, 3),
            (4, 5, 6)
        ]

        @test decompose(Point3{Float64}, Tessellation(s, 8)) ≈ positions

        _faces = TriangleFace[
            (9, 2, 1), (9, 3, 2), (9, 4, 3), (9, 1, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 4, 8), (3, 8, 7), (4, 1, 5), (4, 5, 8),
            (10, 5, 6), (10, 6, 7), (10, 7, 8), (10, 8, 5)]

        @test _faces == decompose(TriangleFace{Int}, Tessellation(s, 8))

        m = triangle_mesh(Tessellation(s, 8))
        @test m === triangle_mesh(m)
        @test GeometryBasics.faces(m) == decompose(GLTriangleFace, _faces)
        @test GeometryBasics.coordinates(m) ≈ positions

        m = normal_mesh(s) # just test that it works without explicit resolution parameter
        @test hasproperty(m, :position)
        @test hasproperty(m, :normal)
        @test length(faces(m)) == length(faces(m.normal))
        @test faces(m) isa AbstractVector{GLTriangleFace}
        @test faces(m.normal) isa AbstractVector{GLTriangleFace}

        ns = GeometryBasics.FaceView(
            Vec{3, Float32}[
                [0.70710677, -0.70710677, 0.0], [0.4082483, 0.4082483, -0.8164966],
                [-0.70710677, 0.70710677, -9.9991995f-17], [-0.4082483, -0.4082483, 0.8164966],
                [-0.57735026, -0.57735026, -0.57735026], [0.57735026, 0.57735026, 0.57735026]
            ], [
                GLTriangleFace(5, 5, 5), GLTriangleFace(5, 5, 5),
                GLTriangleFace(5, 5, 5), GLTriangleFace(5, 5, 5),
                QuadFace{Int64}(1, 2, 2, 1), QuadFace{Int64}(2, 3, 3, 2),
                QuadFace{Int64}(3, 4, 4, 3), QuadFace{Int64}(4, 1, 1, 4),
                GLTriangleFace(6, 6, 6), GLTriangleFace(6, 6, 6),
                GLTriangleFace(6, 6, 6), GLTriangleFace(6, 6, 6)
            ]
        )

        @test ns == decompose_normals(Tessellation(s, 8))

        muv = uv_mesh(s)
        @test !hasproperty(muv, :uv) # not defined yet
    end
end

@testset "HyperRectangles" begin
    @testset "Constructors" begin
        # TODO: Do these actually make sense?
        # Should they not be Rect(NaN..., 0...)?
        @testset "Empty Constructors" begin
            for constructor in [Rect, Rect{2}, Rect2, RectT, Rect2f]
                @test constructor() == Rect{2, Float32}(Inf, Inf, -Inf, -Inf)
            end
            for constructor in [Rect{3}, Rect3, Rect3f]
                @test constructor() == Rect{3, Float32}((Inf, Inf, Inf), (-Inf, -Inf, -Inf))
            end

            for T in [UInt32, Int16, Float64]
                a = typemax(T)
                b = typemin(T)
                for constructor in [Rect{2, T}, Rect2{T}, RectT{T, 2}]
                    @test constructor() == Rect{2, T}(a, a, b, b)
                end
                for constructor in [Rect{3, T}, Rect3{T}, RectT{T, 3}]
                    @test constructor() == Rect{3, T}(Point(a, a, a), Vec(b, b, b))
                end
            end
        end

        @testset "Constructor arg conversions" begin
            function expected_rect(::Type{<: Rect}, arg1, arg2)
                return Rect{min(length(arg1), length(arg2)), promote_type(eltype(arg1), eltype(arg2))}(arg1, arg2)
            end
            function expected_rect(::Type{<: Rect{N}}, arg1, arg2) where {N}
                return Rect{N, promote_type(eltype(arg1), eltype(arg2))}(arg1, arg2)
            end
            function expected_rect(::Type{<: Rect{N, T}}, arg1, arg2) where {N, T}
                return Rect{N, T}(arg1, arg2)
            end

            @testset "2D args -> 2D Rect" begin
                for constructor in [Rect, RectT,   Rect2, Rect{2}, RectT{Int32},
                                    Rect2f, Rect{2, Float16}, Rect2{UInt32}, RectT{Float64, 2}]
                    @testset "$constructor" begin
                        @test constructor(1,2,3,4)           == expected_rect(constructor, Point(1,2), Vec(3,4))
                        @test constructor(1.0,2,3,4)         == expected_rect(constructor, Point(1.0,2), Vec(3,4))
                        @test constructor(Point2f(1,2),3,4)  == expected_rect(constructor, Point2f(1,2), Vec(3,4))
                        @test constructor(Vec2(1,2),3,4.0)   == expected_rect(constructor, Point(1,2), Vec(3,4.0))
                        @test constructor((1,2),Point2(3,4)) == expected_rect(constructor, Point(1,2), Vec(3,4))
                        @test constructor(1.0,2,Vec2(3,4))   == expected_rect(constructor, Point(1,2), Vec(3,4))
                        @test_throws ArgumentError constructor(1,2,3)
                    end
                end
            end

            @testset "3D args -> 3D Rect" begin
                for constructor in [Rect, RectT,   Rect3, Rect{3}, RectT{Float64},
                                    Rect3d, Rect{3, Int16}, Rect3{UInt8}, RectT{Float32, 3}]
                    @testset "$constructor" begin
                        @test constructor(1,2,3,4,5,6)           == expected_rect(constructor, Point(1,2,3), Vec(4,5,6))
                        @test constructor(1,2,3,4,5,6.0)         == expected_rect(constructor, Point(1,2,3), Vec(4,5,6.0))
                        @test constructor(1,2.0,3,Vec3f(4,5,6))  == expected_rect(constructor, Point(1,2,3), Vec3f(4,5,6))
                        @test constructor(Vec3(1,2,3),4,5,6)     == expected_rect(constructor, Point3(1,2,3), Vec(4,5,6))
                        @test constructor((1,2,3),Point3(4,5,6)) == expected_rect(constructor, Point(1,2,3), Vec(4,5,6))
                        @test constructor(Vec3(1,2,3),4,5,6)     == expected_rect(constructor, Point(1,2,3), Vec(4,5,6))
                        @test_throws ArgumentError constructor(1,2,3)
                    end
                end
            end
        end

        @testset "Copy Constructors" begin
            r = Rect2i(0,0,1,1)
            for constructor in [Rect, Rect2f, Rect3f, RectT{Float64}]
                @test constructor(r) == constructor(Point2(0), Vec2(1))
            end
        end

        @testset "Special Constructors" begin
            @test Rect3f((1, 2, 3, Vec(1,2,3))) == Rect3f(1,2,3, Vec(1,2,3))
            @test Rect2(((1, 2), 3, 4)) == Rect2f((1,2), 3, 4)
            @test Rect((1, 2, 3, 4)) == Rect2f(1, 2, 3, 4)
            @test Rect2((x = 1, y = 2), (width = 3, height = 4)) == Rect2f(1, 2, 3, 4)

            for constructor in [Rect, RectT,   Rect2, Rect{2}, RectT{Float32},
                                Rect2f, Rect{2, Float16}, Rect2{Float16}, RectT{Float64, 2}]
                @test centered(constructor) == constructor(Point2d(-0.5), Vec2d(1))
            end
            for constructor in [Rect3, Rect{3}, Rect3d, Rect{3, Float16}, Rect3{Float64}, RectT{Float32, 3}]
                @test centered(constructor) == constructor(Point3d(-0.5), Vec3d(1))
            end
        end

        @testset "Boundingbox-like" begin
            for constructor in [Rect, Rect{2}, Rect{2, Float32}, Rect3f]
                @test constructor(Circle(Point2f(0), 1f0))              == constructor(Point2f(-1, -1), Vec2f(2, 2))
                @test constructor(Rect2f(0, 0, 1, 1))                   == constructor(Point2f( 0,  0), Vec2f(1, 1))
                m = GeometryBasics.mesh(Tessellation(Circle(Point2f(0), 1f0), 5))
                @test constructor(m)                                     ≈ constructor(Point2f(-1, -1), Vec2f(2, 2))
            end
            for constructor in [Rect, Rect{3}, Rect{3, Float32}]
                @test constructor(Sphere(Point3f(0), 1f0))                          == Rect3f(-1, -1, -1, 2, 2, 2)
                @test constructor(Rect3f(0, 0, 0, 1, 1, 1))                         == Rect3f(0, 0, 0, 1, 1, 1)
                @test constructor(Cylinder(Point3f(0, 0, -1), Point3f(0,0,1), 1f0))  ≈ Rect3f(-1, -1, -1, 2, 2, 2) atol = 0.05
                @test constructor(Pyramid(Point3f(0, 0, -1), 2f0, 2f0))             == Rect3f(-1, -1, -1, 2, 2, 2)
                m = GeometryBasics.mesh(Tessellation(Sphere(Point3f(0), 1f0), 5))
                @test constructor(m)                                                 ≈ Rect3f(-1, -1, -1, 2, 2, 2)
            end
        end
    end

    # TODO: consider deprecating this
    @testset "Integer rounding" begin
        @test Rect2i(0.3, 0.4, 0.6, 0.5) == Rect2i(0,0,1,0)
        @test Rect2{UInt32}(0.8, 0.1, 1.3, 1.9) == Rect2i(1,0,1,2)
        @test Rect3i(0.3, 0.6, 0.9, 1.2, 1.5, 1.8) == Rect3i(0,1,1,1,2,2)
    end
    @testset "Vector inputs" begin
        @test Rect([1,2], [2,3]) == Rect2i(1,2,2,3)
        @test Rect2f([0,1], [2,3]) == Rect2f(0,1,2,3)
        @test RectT{Float64}([0,1,0], [2,2,2]) == Rect3d(0,1,0, 2,2,2)
        @test Rect3([0,1,0], [2,2,2]) == Rect3i(0,1,0, 2,2,2)
    end

    # TODO: These don't really make sense...
    r = Rect2f()
    @test origin(r) == Vec(Inf, Inf)
    @test minimum(r) == Vec(Inf, Inf)
    @test isnan(maximum(r))
    @test width(r) == -Inf
    @test height(r) == -Inf
    @test widths(r) == Vec(-Inf, -Inf)
    @test area(r) == Inf
    @test volume(r) == Inf
    # TODO: broken? returns NaN widths
    # @test union(r, Rect2f(1,1,2,2)) == Rect2f(1,1,2,2)
    # @test union(Rect2f(1,1,2,2), r) == Rect2f(1,1,2,2)
    @test update(r, Vec2f(1,1)) == Rect2f(1,1,0,0)

    a = Rect(Vec(0, 1), Vec(2, 3))
    pt_expa = Point{2,Int}[(0, 1), (2, 1), (2, 4), (0, 4)]
    @test decompose(Point{2,Int}, a) == pt_expa
    mesh = normal_mesh(a)
    @test decompose(Point2f, mesh) == pt_expa

    @test origin(a) == Vec(0,1)
    @test minimum(a) == Vec(0,1)
    @test maximum(a) == Vec(2,4)
    @test width(a) == 2
    @test height(a) == 3
    @test widths(a) == Vec(2,3)
    @test area(a) == 2*3
    @test volume(a) == 2*3
    @test union(a, Rect2f(1,1,2,2)) == Rect2f(0,1,3,3)
    @test union(Rect2f(1,1,2,2), a) == Rect2f(0,1,3,3)
    @test update(a, Vec2f(0,0)) == Rect2f(0,0,2,4)

    b = Rect(Vec(1, 2, 3), Vec(4, 5, 6))
    pt_expb = Point{3, Int64}[[1, 2, 3], [1, 2, 9], [1, 7, 3], [1, 7, 9],
                              [5, 2, 3], [5, 2, 9], [5, 7, 3], [5, 7, 9]]
    @test decompose(Point{3,Int}, b) == pt_expb

    @test origin(b) == Vec(1,2,3)
    @test minimum(b) == Vec(1,2,3)
    @test maximum(b) == Vec(5,7,9)
    @test width(b) == 4
    @test height(b) == 5
    @test widths(b) == Vec(4,5,6)
    @test_throws MethodError area(b)
    @test volume(b) == 4*5*6
    @test union(b, Rect3f(1,1,1,2,2,2)) == Rect3f(1,1,1, 4,6,8)
    @test union(Rect3f(1,1,1,2,2,2), b) == Rect3f(1,1,1, 4,6,8)
    @test update(b, Vec3f(0)) == Rect3f(0,0,0,5,7,9)


    mesh = normal_mesh(b)
    @test faces(mesh) == GLTriangleFace[
        (1, 2, 4), (1, 4, 3), (7, 8, 6), (7, 6, 5), (5, 6, 2), (5, 2, 1),
        (3, 4, 8), (3, 8, 7), (1, 3, 7), (1, 7, 5), (6, 8, 4), (6, 4, 2)]
    @test normals(mesh) == GeometryBasics.FaceView(
        Vec{3, Float32}[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0]],
        GLTriangleFace[(1, 1, 1), (1, 1, 1), (2, 2, 2), (2, 2, 2), (3, 3, 3), (3, 3, 3), (4, 4, 4), (4, 4, 4), (5, 5, 5), (5, 5, 5), (6, 6, 6), (6, 6, 6)]
    )
    @test coordinates(mesh) == Point{3, Float32}[
        [1, 2, 3], [1, 2, 9], [1, 7, 3], [1, 7, 9],
        [5, 2, 3], [5, 2, 9], [5, 7, 3], [5, 7, 9]]

    @test isempty(Rect{3,Float32}())
end

@testset "Pyramids" begin
    p = Pyramid(Point3f(0), 1f0, 0.2f0)
    @test coordinates(p) == Point3f[[0.0, 0.0, 1.0], [0.1, -0.1, 0.0], [0.1, 0.1, 0.0], [-0.1, 0.1, 0.0], [-0.1, -0.1, 0.0]]
    @test faces(p) == [GLTriangleFace(1, 2, 3), GLTriangleFace(1, 3, 4), GLTriangleFace(1, 4, 5), GLTriangleFace(1, 5, 2), QuadFace{Int64}(2, 3, 4, 5)]
    ns = normals(p)
    @test faces(ns) == [GLTriangleFace(1), GLTriangleFace(2), GLTriangleFace(3), GLTriangleFace(4), QuadFace{Int64}(5)]
    @test values(ns) ≈ Vec3f[[0.9805807, 0.0, 0.19611615], [0.0, 0.9805807, 0.19611615], [-0.9805807, 0.0, 0.19611615], [0.0, -0.9805807, 0.19611615], [0.0, 0.0, -1.0]]
end

@testset "Type Utils" begin
    @test GeometryBasics.Face(TriangleFace, Int32) == TriangleFace{Int32}
    @test GeometryBasics.Face(GLTriangleFace, Int32) == TriangleFace{GLIndex}
end

@testset "Faces" begin
    @test convert_simplex(GLTriangleFace, QuadFace{Int}(1, 2, 3, 4)) ==
          (GLTriangleFace(1, 2, 3), GLTriangleFace(1, 3, 4))
    @test convert_simplex(NgonFace{3,ZeroIndex{Int}}, QuadFace{ZeroIndex{Int}}(1, 2, 3, 4)) ==
          (NgonFace{3,ZeroIndex{Int}}(1, 2, 3), NgonFace{3,ZeroIndex{Int}}(1, 3, 4))
    @test convert_simplex(NgonFace{3,OffsetInteger{3,Int}},
                          NgonFace{4,OffsetInteger{2,Int}}(1, 2, 3, 4)) ==
          (NgonFace{3,OffsetInteger{3,Int}}(1, 2, 3), NgonFace{3,OffsetInteger{3,Int}}(1, 3, 4))
    @test convert_simplex(LineFace{Int}, QuadFace{Int}(1, 2, 3, 4)) ==
          (LineFace{Int}(1, 2), LineFace{Int}(2, 3), LineFace{Int}(3, 4),
           LineFace{Int}(4, 1))

    @testset "NgonFace ambiguity" begin
        face = NgonFace((1, 2))
        @test convert_simplex(NgonFace{2,UInt32}, face) === (NgonFace{2,UInt32}((1, 2)),)
        @test convert_simplex(typeof(face), face) === (face,)
        face = NgonFace((1,))
        @test convert_simplex(NgonFace{1,UInt32}, face) === (NgonFace{1,UInt32}((1,)),)
        @test convert_simplex(typeof(face), face) === (face,)
    end

    ps = rand(Point2f, 10)
    f = GLTriangleFace(1, 2, 3)
    @test ps[f] == Triangle(ps[[1,2,3]]...)
    data = [string(i) for i in 1:10]
    f = QuadFace(3, 4, 7, 8)
    @test data[f] == ("3", "4", "7", "8")

    @test GeometryBasics.cyclic_hash(f) != GeometryBasics.cyclic_hash(QuadFace(1,2,3,4))
    @test GeometryBasics.cyclic_hash(f) == GeometryBasics.cyclic_hash(QuadFace(3,4,7,8))
    # cyclic permutation does not change the face
    @test GeometryBasics.cyclic_hash(f) == GeometryBasics.cyclic_hash(QuadFace(7,8,3,4))
    @test GeometryBasics.cyclic_hash(GLTriangleFace(1,2,3)) == GeometryBasics.cyclic_hash(GLTriangleFace(1,2,3))
    @test GeometryBasics.cyclic_hash(GLTriangleFace(1,2,3)) == GeometryBasics.cyclic_hash(GLTriangleFace(2,3,1))
    @test GeometryBasics.cyclic_hash(GLTriangleFace(1,2,3)) == GeometryBasics.cyclic_hash(GLTriangleFace(3,1,2))

    # repeat with cyclic_equal
    @test !GeometryBasics.cyclic_equal(f, QuadFace(1,2,3,4))
    @test GeometryBasics.cyclic_equal(f, QuadFace(3,4,7,8))
    @test GeometryBasics.cyclic_equal(f, QuadFace(7,8,3,4))
    @test GeometryBasics.cyclic_equal(GLTriangleFace(1,2,3), GLTriangleFace(1,2,3))
    @test GeometryBasics.cyclic_equal(GLTriangleFace(1,2,3), GLTriangleFace(2,3,1))
    @test GeometryBasics.cyclic_equal(GLTriangleFace(1,2,3), GLTriangleFace(3,1,2))
end

@testset "FaceView" begin
    ps = rand(Point2f, 5)
    fs = GLTriangleFace[(1,2,3), (2,3,4), (5,5,5)]
    fv = FaceView(ps, fs)
    @test faces(fv) == fs
    @test values(fv) == ps
    @test fv[fs[1]] == ps[fs[1]]
    @test !isempty(fv)
    @test fv == FaceView(ps, fs)
    @test length(fv) == 5
end

@testset "Normals" begin
    # per face normals
    r = Rect3f(Point3f(0), Vec3f(1))
    ns = face_normals(coordinates(r), faces(r))
    ux = unit(Vec3f, 1); uy = unit(Vec3f, 2); uz = unit(Vec3f, 3)
    @test ns == normals(r)
    @test values(ns) == [-ux, ux, -uy, uy, -uz, uz]

    # typing
    ux = unit(Vec3d, 1); uy = unit(Vec3d, 2); uz = unit(Vec3d, 3)
    ns = face_normals(decompose(Point3d, r), faces(r))
    @test ns isa FaceView{Vec3d}
    @test values(ns) == [-ux, ux, -uy, uy, -uz, uz]

    # Mixed
    c = Cylinder(Point3f(0), Point3f(0,0,1), 0.5f0)
    ns = normals(c)
    # caps without mantle
    f_ns = face_normals(coordinates(c), filter!(f -> f isa TriangleFace, faces(c)))
    @test all(n -> n == values(ns)[end-1], values(f_ns)[1:15])
    @test all(n -> n == values(ns)[end], values(f_ns)[16:end])
    # Mantle without caps
    v_ns = normals(coordinates(c), filter!(f -> f isa QuadFace, faces(c)))[1:end-2]
    @test values(ns)[1:15] ≈ v_ns[1:15]
    @test values(ns)[1:15] ≈ v_ns[16:30] # repeated via FaceView in ns
end

@testset "HyperSphere" begin
    sphere = Sphere{Float32}(Point3f(0), 1.0f0)

    points = decompose(Point, Tessellation(sphere, 3))
    point_target = Point{3,Float32}[[0.0, 0.0, 1.0], [1.0, 0.0, 6.12323e-17],
                                    [1.22465e-16, 0.0, -1.0], [-0.0, 0.0, 1.0],
                                    [-1.0, 1.22465e-16, 6.12323e-17],
                                    [-1.22465e-16, 1.49976e-32, -1.0], [0.0, -0.0, 1.0],
                                    [1.0, -2.44929e-16, 6.12323e-17],
                                    [1.22465e-16, -2.99952e-32, -1.0]]
    @test points ≈ point_target

    f = decompose(TriangleFace{Int}, Tessellation(sphere, 3))
    face_target = TriangleFace{Int}[[1, 2, 5], [1, 5, 4], [2, 3, 6], [2, 6, 5], [4, 5, 8],
                                    [4, 8, 7], [5, 6, 9], [5, 9, 8]]
    @test f == face_target

    uv = decompose_uv(Tesselation(sphere, 3))
    uv_target = Vec{2, Float32}[[0.0, 1.0], [0.0, 0.5], [0.0, 0.0], [0.5, 1.0], [0.5, 0.5],
                                [0.5, 0.0], [1.0, 1.0], [1.0, 0.5], [1.0, 0.0]]
    @test uv == uv_target

    @test minimum(sphere) == Point3f(-1)
    @test maximum(sphere) == Point3f(1)
    @test origin(sphere) == Point3f(0)
    @test widths(sphere) == Vec3f(2)
    @test radius(sphere) == 1f0
    @test !(Point3f(1) in sphere)
    @test Point3f(0.5) in sphere
    @test centered(HyperSphere) == Sphere(Point3f(0), 0.5f0)
    @test centered(Sphere) == Sphere(Point3f(0), 0.5f0)
    @test centered(Sphere{Float64}) == Sphere(Point3(0.0), 0.5)

    circle = Circle(Point2f(0), 1.0f0)
    points = decompose(Point2f, Tessellation(circle, 20))
    @test length(points) == 20
    tess_circle = Tessellation(circle, 32)
    mesh = triangle_mesh(tess_circle)
    @test decompose(Point2f, mesh) ≈ decompose(Point2f, tess_circle)

    @test minimum(circle) == Point2f(-1)
    @test maximum(circle) == Point2f(1)
    @test origin(circle) == Point2f(0)
    @test widths(circle) == Vec2f(2)
    @test radius(circle) == 1f0
    @test !(Point2f(-1) in circle)
    @test Point2f(-0.5) in circle
    @test centered(Circle) == Circle(Point2f(0), 0.5f0)
    @test centered(Circle{Float64}) == Circle(Point2(0.0), 0.5)
end

@testset "Rectangles" begin
    rect = Rect2f(0, 7, 20, 3)
    @test (rect + 4) == Rect2f(4, 11, 20, 3)
    @test (rect + Vec(2, -2)) == Rect2f(2, 5, 20, 3)

    @test (rect - 4) == Rect2f(-4, 3, 20, 3)
    @test (rect - Vec(2, -2)) == Rect2f(-2, 9, 20, 3)

    base = Vec3f(1, 2, 3)
    wxyz = Vec3f(-2, 4, 2)
    rect = Rect3f(base, wxyz)
    @test (rect + 4) == Rect3f(base .+ 4, wxyz)
    @test (rect + Vec(2, -2, 3)) == Rect3f(base .+ Vec(2, -2, 3), wxyz)

    @test (rect - 4) == Rect3f(base .- 4, wxyz)
    @test (rect - Vec(2, -2, 7)) == Rect3f(base .- Vec(2, -2, 7), wxyz)

    rect = Rect2f(0, 7, 20, 3)
    @test (rect * 4) == Rect2f(0, 7 * 4, 20 * 4, 3 * 4)
    @test (rect * Vec(2, -2)) == Rect2f(0, -7 * 2, 20 * 2, -3 * 2)

    base = Vec3f(1, 2, 3)
    wxyz = Vec3f(-2, 4, 2)
    rect = Rect3f(base, wxyz)
    @test (rect * 4) == Rect3f(base .* 4, wxyz .* 4)
    @test (rect * Vec(2, -2, 3)) == Rect3f(base .* Vec(2, -2, 3), wxyz .* Vec(2, -2, 3))

    rect1 = Rect(Vec(0.0, 0.0), Vec(1.0, 2.0))
    rect2 = Rect(0.0, 0.0, 1.0, 2.0)
    @test rect1 isa GeometryBasics.HyperRectangle{2,Float64}
    @test rect1 == rect2

    split1, split2 = GeometryBasics.split(rect1, 2, 1)
    @test widths(split1) == widths(split2)
    @test origin(split1) == Vec(0, 0)
    @test origin(split2) == Vec(0, 1)
    @test in(split1, rect1) && in(split2, rect1)
    @test !(in(rect1, split1) || in(rect1, split2))

    rect1 = Rect(Vec(0.0, 0.0, -1.0), Vec(1.0, 2.0, 1.0))
    split1, split2 = GeometryBasics.split(rect1, 1, 0.75)
    @test widths(split1) == Vec(0.75, 2, 1)
    @test widths(split2) == Vec(0.25, 2, 1)
    @test origin(split1) == Vec(0, 0, -1)
    @test origin(split2) == Vec(0.75, 0, -1)
    @test in(split1, rect1) && in(split2, rect1)
    @test !(in(rect1, split1) || in(rect1, split2))

    prim = Rect(0.0, 0.0, 1.0, 1.0)
    @test length(prim) == 2

    @test width(prim) == 1.0
    @test height(prim) == 1.0

    b1 = Rect2(0.0, 0.0, 2.0, 2.0)
    b2 = Rect2(0, 0, 2, 2)
    @test isequal(b1, b2)

    pt = Point(1.0, 1.0)
    b1 = Rect(0.0, 0.0, 1.0, 1.0)
    @test in(pt, b1)

    rect = Rect(0.0, 0.0, 1.0, 1.0)
    @test GeometryBasics.positive_widths(rect) isa GeometryBasics.HyperRectangle{2,Float64}

    h1 = Rect(0.0, 0.0, 1.0, 1.0)
    h2 = Rect(1.0, 1.0, 2.0, 2.0)
    @test union(h1, h2) isa GeometryBasics.HyperRectangle{2,Float64}
    # @test GeometryBasics.diff(h1, h2) == h1
    @test GeometryBasics.intersect(h1, h2) isa GeometryBasics.HyperRectangle{2,Float64}

    b = Rect(0.0, 0.0, 1.0, 1.0)
    v = Vec(1, 2)
    @test update(b, v) isa GeometryBasics.HyperRectangle{2,Float64}
    v = Vec(1.0, 2.0)
    @test update(b, v) isa GeometryBasics.HyperRectangle{2,Float64}

    @testset "euclidean distances" begin
        p = Vec(5.0, 4.0)
        rect = Rect(0.0, 0.0, 1.0, 1.0)
        @test min_dist_dim(rect, p, 1) == 4.0
        @test min_dist_dim(rect, p, 2) == 3.0
        @test max_dist_dim(rect, p, 1) == 5.0
        @test max_dist_dim(rect, p, 2) == 4.0
        @test minmax_dist_dim(rect, p, 1) == (4.0, 5.0)

        rect1 = Rect(0.0, 0.0, 1.0, 1.0)
        rect2 = Rect(3.0, 1.0, 4.0, 2.0)
        @test min_dist_dim(rect1, rect2, 1) == 2.0
        @test min_dist_dim(rect1, rect2, 2) == 0.0
        @test max_dist_dim(rect1, rect2, 1) == 7.0
        @test max_dist_dim(rect1, rect2, 2) == 3.0
        @test minmax_dist_dim(rect1, rect2, 1) == (2.0, 7.0)

        r = Rect2f(-1, -1, 2, 3)
        p = Point2f(1, 2) + Point2f(3, 4)
        @test min_euclidean(r, p) == 5f0
        @test max_euclidean(r, p) ≈ sqrt(5*5 + 7*7)

        r2 = Rect2f(0, 0, 2, 3)
        @test min_euclidean(r, r2) == 0f0
        @test max_euclidean(r, r2) == 5f0
        @test minmax_euclidean(r, r2) == (0f0, 5f0)
    end

    @test !before(rect1, rect2)
    rect1 = Rect(0.0, 0.0, 1.0, 1.0)
    rect2 = Rect(3.0, 2.0, 4.0, 2.0)
    @test before(rect1, rect2)

    @test !meets(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 4.0, 2.0)
    @test meets(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !overlaps(rect1, rect2)
    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(1.5, 1.5, 2.0, 2.0)
    @test overlaps(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !GeometryBasics.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 1.5, 1.5)
    @test !GeometryBasics.starts(rect1, rect2)
    rect2 = Rect(1.0, 1.0, 3.0, 3.0)
    @test GeometryBasics.starts(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 4.0, 4.0)
    @test during(rect1, rect2)
    rect1 = Rect(0.0, 0.0, 2.0, 3.0)
    rect2 = Rect(1.0, 1.0, 4.0, 2.0)
    @test !during(rect1, rect2)

    rect1 = Rect(1.0, 1.0, 2.0, 2.0)
    rect2 = Rect(0.0, 0.0, 4.0, 4.0)
    @test !finishes(rect1, rect2)
    rect1 = Rect(1.0, 0.0, 1.0, 1.0)
    rect2 = Rect(0.0, 0.0, 2.0, 1.0)
    @test !finishes(rect1, rect2)
    rect1 = Rect(1.0, 1.0, 1.0, 2.0)
    rect2 = Rect(0.0, 0.0, 2.0, 3.0)
    @test finishes(rect1, rect2)

    rect1 = @inferred Rect(1, 2, 3, 4, 5, 6, 7, 8)
    rect2 = Rect(Vec(1, 2, 3, 4), Vec(5, 6, 7, 8))
    @test rect1 == rect2

    @testset "Matrix Multiplications" begin
        r = Rect2f(-1, -2, 4, 3)

        # TODO: this seems quite dangerous: We pad points with ones which makes
        #       sense for translations if we go to D+1, but is nonsense if we
        #       go higher dimensions than that.
        M = rand(Mat4f)
        ps = Point2f[M * Point(p..., 1, 1) for p in coordinates(r)]
        @test Rect2f(ps) ≈ M * r

        M = Mat2f(0.5, -0.3, 0.7, 1.5)
        ps = Point2f[M * p for p in coordinates(r)]
        @test Rect2f(ps) ≈ M * r

        r = Rect3f(-1, -2, -3, 2, 4, 1)
        M = rand(Mat4f)
        ps = Point3f[M * Point(p..., 1) for p in coordinates(r)]
        @test Rect3f(ps) ≈ M * r
    end

    # TODO: this is effectively 0-indexed... should it be?
    M = reshape(collect(11:100), 10, 9)[1:9, :]
    r = Rect2i(2, 4, 2, 4)
    @test M[r] == [53 63 73 83; 54 64 74 84]

end

@testset "LineStrings" begin
    ps1 = rand(Point2f, 10)
    ls1 = LineString(ps1)
    _ls1 = LineString(ps1)
    @test coordinates(ls1) == ps1
    @test length(ls1) == 10
    @test ls1 == _ls1

    ls2 = LineString(rand(Point2f, 6))
    ls3 = LineString(rand(Point2f, 4))
    mls = MultiLineString([ls1, ls2, ls3])
    @test mls.linestrings == [ls1, ls2, ls3]
    @test mls[1] == ls1
    @test mls[2] == ls2
    @test mls[3] == ls3
    @test size(mls) == (3, ) # TODO: Does this make sense?
    @test length(mls) == 3
    @test MultiLineString(OffsetArray([ls1, ls2, ls3], 0)) == mls
end

@testset "MultiPoint" begin
    ps1 = rand(Point2f, 10)
    mp = MultiPoint(ps1)
    @test all(getindex.(Ref(mp), 1:10) .== ps1)
    @test size(mp) == (10, ) # TODO: Does this make sense?
    @test length(mp) == 10
end