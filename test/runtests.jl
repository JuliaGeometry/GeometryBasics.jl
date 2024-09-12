using Test, Random, OffsetArrays
using GeometryBasics
using LinearAlgebra
using GeometryBasics: MetaMesh
using GeoInterface
using GeoJSON
using Extents

@testset "GeometryBasics" begin
@testset "algorithms" begin
    cube = Rect(Vec3f(-0.5), Vec3f(1))
    cube_faces = decompose(TriangleFace{Int}, faces(cube))
    cube_vertices = decompose(Point{3,Float32}, cube)
    @test area(cube_vertices, cube_faces) == 6
    mesh = Mesh(cube_vertices, cube_faces)
    @test GeometryBasics.volume(mesh) ≈ 1
    @test GeometryBasics.volume(cube) ≈ 1

    rect = Rect(1, 2, 7.5, 2.0)
    @test GeometryBasics.area(rect) ≈ 15

    points_cwise = Point2f[(0,0), (0,1), (1,1)]
    points_ccwise = Point2f[(0,0), (1,0), (1,1)]
    @test area(points_cwise) ≈ -0.5
    @test area(points_ccwise) ≈ 0.5
    @test area(OffsetArray(points_cwise, -2)) ≈ -0.5

    points3d = Point3f[(0,0,0), (0,0,1), (0,1,1)]
    @test area(OffsetArray(points3d, -2)) ≈ 0.5
end

@testset "embedding metadata" begin
    @testset "Meshes" begin
        @testset "per vertex attributes" begin
            points = rand(Point{3, Float64}, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(Vec{3, Float64}, 8)
            stress = LinRange(0, 1, 8)
            mesh = Mesh(points, tfaces; normal = normals, stress = stress)

            @test hasproperty(mesh, :stress)
            @test hasproperty(mesh, :normal)
            @test mesh.stress === stress
            @test mesh.normal === normals
            @test mesh.position === points
            @test GeometryBasics.faces(mesh) === tfaces
            # TODO: Is the order variable with Dict?
            @test propertynames(mesh) == (:vertex_attributes, :connectivity, :views, :normal, :position, :stress)
        end
    end

    @testset "Mesh with metadata" begin
       m = triangle_mesh(Sphere(Point3f(0), 1))
       m_meta = MetaMesh(m; boundingbox=Rect(1.0, 1.0, 2.0, 2.0))
       @test m_meta[:boundingbox] === Rect(1.0, 1.0, 2.0, 2.0)
       @test collect(keys(m_meta)) == [:boundingbox,]
    end
end

@testset "view" begin
    @testset "TupleView" begin
        x = [1, 2, 3, 4, 5, 6]
        y = TupleView{2, 1}(x)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]

        y = TupleView{2}(x)
        @test y == [(1, 2), (3, 4), (5, 6)]

        y = TupleView{2, 3}(x)
        @test y == [(1, 2), (4, 5)]

        y = TupleView{3, 1}(x)
        @test y == [(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6)]

        y = TupleView{2, 1}(x, connect = true)
        @test y == [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]

    end

    @testset "connected views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        x = connect(numbers, Point{2})

        @test x == Point[(1, 2), (3, 4), (5, 6)]

        line = connect(x, Line, 1)
        @test line == [Line(Point(1, 2), Point(3, 4)), Line(Point(3, 4), Point(5, 6))]

        triangles = connect(x, Triangle)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
        x = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        tetrahedra = connect(x, NSimplex{4})
        @test tetrahedra == [Tetrahedron(x[1], x[2], x[3], x[4])]

        @testset "matrix non-copy point views" begin
            # point in row
            points = [1 2; 1 4; 66 77]
            comparison = [Point(1, 2), Point(1, 4), Point(66, 77)]
            @test connect(points, Point{2}) == comparison
            # point in column
            points = [1 1 66; 2 4 77]
            # huh, reinterpret array doesn't seem to like `==`
            @test all(((a,b),)-> a==b, zip(connect(points, Point{2}), comparison))
        end
    end

    @testset "reinterpret" begin
        numbers = collect(reshape(1:6, 2, 3))
        points = reinterpret(Point{2, Int}, numbers)
        @test points[1] === Point(1, 2)
        @test points[2] === Point(3, 4)
        numbers[4] = 0
        @test points[2] === Point(3, 0)
    end
end

@testset "constructors" begin
    @testset "Mesh" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})

        mesh = Mesh(points, [1,2,3])
        @test faces(mesh) == [TriangleFace(1, 2, 3)]

        x = Point{3}(1.0)
        mesh = Mesh([x], [TriangleFace(1, 1, 1)])
        @test coordinates(mesh) == [x]
        @test faces(mesh) == [TriangleFace(1, 1, 1)]

        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        f = connect([1, 2, 3, 4], SimplexFace{4})
        mesh = Mesh(points, f)
        @test collect(mesh) == [Tetrahedron(points...)]

        points = rand(Point3f, 8)
        tfaces = [GLTriangleFace(1, 2, 3), GLTriangleFace(5, 6, 7)]
        ns = rand(Vec3f, 8)
        uv = rand(Vec2f, 8)
        mesh = Mesh(points, tfaces)
        meshuv = MetaMesh(points, tfaces; uv=uv)
        meshuvnormal = MetaMesh(points, tfaces; normal=ns, uv=uv)
        t = Tesselation(Rect2f(0, 0, 2, 2), (30, 30))

        m = GeometryBasics.mesh(t; pointtype=Point3f, facetype=QuadFace{Int})
        m2 = GeometryBasics.mesh(m, facetype=QuadFace{GLIndex})
        @test GeometryBasics.faces(m2) isa Vector{QuadFace{GLIndex}}
        @test GeometryBasics.coordinates(m2) isa Vector{Point3f}
    end
end

@testset "decompose/triangulation" begin
    @test isempty(decompose(Vec3f, []))
    @test decompose(Vec3f, []) isa Vector{Vec3f}
    primitive = Sphere(Point3f(0), 1)
    @test ndims(primitive) === 3
    mesh = triangle_mesh(primitive)
    @test decompose(Point, mesh) isa Vector{Point3f}
    @test decompose(Point, primitive) isa Vector{Point3f}

    primitive = Rect2(0, 0, 1, 1)
    mesh = triangle_mesh(primitive)

    @test decompose(Point, mesh) isa Vector{Point2f}
    @test decompose(Point, primitive) isa Vector{Point2{Int}}

    primitive = Rect3(0, 0, 0, 1, 1, 1)
    triangle_mesh(primitive)

    primitive = Sphere(Point3f(0), 1)
    m_normal = normal_mesh(primitive)
    @test GeometryBasics.normals(m_normal) isa Vector{Vec3f}
    primitive = Rect2(0, 0, 1, 1)
    m_normal = normal_mesh(primitive)
    @test GeometryBasics.normals(m_normal) isa Vector{Vec3f}
    primitive = Rect3(0, 0, 0, 1, 1, 1)
    m_normal = normal_mesh(primitive)
    @test GeometryBasics.normals(m_normal) isa GeometryBasics.FaceView{Vec3f, Vector{Vec3f}, Vector{GLTriangleFace}}

    points = decompose(Point2f, Circle(Point2f(0), 1))
    tmesh = triangle_mesh(points)
    @test GeometryBasics.normals(tmesh) == nothing

    m = GeometryBasics.mesh(Sphere(Point3f(0), 1))
    @test GeometryBasics.normals(m) == nothing
    m_normals = decompose_normals(m)
    @test m_normals isa Vector{Vec3f}

    @test texturecoordinates(m) == nothing
    r2 = Rect2(0.0, 0.0, 1.0, 1.0)
    @test collect(texturecoordinates(r2)) == Point2f[(0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]
    r3 = Rect3(0.0, 0.0, 1.0, 1.0, 2.0, 2.0)
    @test first(texturecoordinates(r3)) == Vec3(0, 0, 0)
    uv = decompose_uv(m)
    @test_broken false # Rect(Point.(uv)) == Rect(0, 0, 1, 1) # decompose_uv must now produces 2D uvs
    uvw = GeometryBasics.decompose_uvw(m)
    @test Rect(Point.(uvw)) == Rect(Point3f(0), Vec3f(1))

    points = decompose(Point2f, Circle(Point2f(0), 1))
    m = GeometryBasics.mesh(points)
    @test coordinates(m) === points
end

@testset "convert mesh + meta" begin
    m = uv_normal_mesh(Circle(Point2f(0), 1f0))
    # For 2d primitives normal is just the upvector
    @test m.normal == [Vec3f(0, 0, 1) for p in coordinates(m)]
end

@testset "mesh conversion" begin
    s = Sphere(Point3(0.0), 1.0)
    m = GeometryBasics.mesh(s)
    @test m isa Mesh{3, Float64}
    @test coordinates(m) isa Vector{Point{3, Float64}}
    @test GeometryBasics.faces(m) isa Vector{GLTriangleFace}
    # Check, that decompose isn't making a copy for matching eltype
    @test coordinates(m) === decompose(Point{3, Float64}, m)

    tmesh = triangle_mesh(m)
    @test coordinates(tmesh) === decompose(Point3f, tmesh)

    nmesh = normal_mesh(m)
    @test coordinates(nmesh) === decompose(Point3f, nmesh)
    @test GeometryBasics.normals(nmesh) === decompose_normals(nmesh)

    m = GeometryBasics.mesh(s, pointtype=Point3f)
    @test m isa Mesh{3, Float32}
    @test coordinates(m) isa Vector{Point3f}
    @test GeometryBasics.faces(m) isa Vector{GLTriangleFace}
end

@testset "lines intersects" begin
    a = Line(Point(0.0, 0.0), Point(4.0, 1.0))
    b = Line(Point(0.0, 0.25), Point(3.0, 0.25))
    c = Line(Point(0.0, 0.25), Point(0.5, 0.25))
    d = Line(Point(0.0, 0.0), Point(0.0, 4.0))
    e = Line(Point(1.0, 0.0), Point(0.0, 4.0))
    f = Line(Point(5.0, 0.0), Point(6.0, 0.0))

    @test intersects(a, b) === (true, Point(1.0, 0.25))
    @test intersects(a, c) === (false, Point(0.0, 0.0))
    @test intersects(d, d) === (false, Point(0.0, 0.0))
    found, point = intersects(d, e)
    @test found && point ≈ Point(0.0, 4.0)
    @test intersects(a, f) === (false, Point(0.0, 0.0))

    # issue #168
    # If these tests fail then you can increase the tolerance on the checks so
    # long as you know what you're doing :)
    line_helper(a, b, c, d) = Line(Point(a, b), Point(c, d))
    b, loc = intersects(line_helper(-3.1, 15.588457268119894, 3.1, 15.588457268119894),
                        line_helper(2.0866025403784354, 17.37050807568877, -4.0866025403784505, 13.806406460551015))
    @test b
    @test loc ≈ Point(-1.0000000000000058, 15.588457268119894)

    b, loc = intersects(line_helper(5743.933982822018, 150.0, 5885.355339059327, -50.0),
                        line_helper(5760.0, 100.0, 5760.0, 140.0))
    @test b
    @test loc ≈ Point(5760.0, 127.27922061357884)
end

@testset "Offsetintegers" begin
    x = 1
    @test GeometryBasics.raw(x) isa Int64
    @test GeometryBasics.value(x) == x

    x = ZeroIndex(1)
    @test eltype(x) == Int64

    x = OffsetInteger{0}(1)
    @test typeof(x) == OffsetInteger{0,Int64}

    x1 = OffsetInteger{0}(2)
    x2 = 1
    @test Base.to_index(x1) == 2
    @test -(x1) == OffsetInteger{0,Int64}(-2)
    @test abs(x1) == OffsetInteger{0,Int64}(2)
    @test +(x, x1) == OffsetInteger{0,Int64}(3)
    @test *(x, x1) == OffsetInteger{0,Int64}(2)
    @test -(x, x1) == OffsetInteger{0,Int64}(-1)
    #test for /
    @test div(x, x1) == OffsetInteger{0,Int64}(0)
    @test !==(x, x1)
    @test !>=(x, x1)
    @test <=(x, x1)
    @test !>(x, x1)
    @test <(x, x1)
end

@testset "Tests from GeometryTypes" begin
    include("geometrytypes.jl")
end

@testset "Point & Vec type" begin
    include("fixed_arrays.jl")
end

@testset "Some mesh issues" begin
    include("meshes.jl")
end

@testset "GeoInterface" begin
    include("geointerface.jl")
end

using Aqua
# Aqua tests
# Intervals brings a bunch of ambiquities unfortunately
# seems like we also run into https://github.com/JuliaTesting/Aqua.jl/issues/86
Aqua.test_all(GeometryBasics; ambiguities=false, unbound_args=false)

end  # testset "GeometryBasics"
