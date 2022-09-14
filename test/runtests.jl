using Test, Random, StructArrays, Tables, StaticArrays, OffsetArrays
using GeometryBasicsCore, GeometryBasics
using GeometryBasics: attributes
using LinearAlgebra
using GeoInterface
using GeoJSON

@testset "GeometryBasics" begin

@testset "algorithms" begin
    cube = Rect(Vec3f(-0.5), Vec3f(1))
    cube_faces = decompose(TriangleFace{Int}, faces(cube))
    cube_vertices = decompose(Point{3,Float32}, cube)
    @test area(cube_vertices, cube_faces) == 6
    mesh = Mesh(cube_vertices, cube_faces)
    @test volume(mesh) ≈ 1
    @test volume(cube) ≈ 1

    rect = Rect(1, 2, 7.5, 2.0)
    @test area(rect) ≈ 15

    points_cwise = Point2f[(0,0), (0,1), (1,1)]
    points_ccwise = Point2f[(0,0), (1,0), (1,1)]
    @test area(points_cwise) ≈ -0.5
    @test area(points_ccwise) ≈ 0.5
    @test area(OffsetArray(points_cwise, -2)) ≈ -0.5

    points3d = Point3f[(0,0,0), (0,0,1), (0,1,1)]
    @test area(OffsetArray(points3d, -2)) ≈ 0.5

    pm2d = [PointMeta(0.0, 0.0, a=:d), PointMeta(0.0, 1.0, a=:e), PointMeta(1.0, 0.0, a=:f)]
    @test area(pm2d) ≈ -0.5

    pm3d = [PointMeta(0.0, 0.0, 0.0, a=:d), PointMeta(0.0, 1.0, 0.0, a=:e), PointMeta(1.0, 0.0, 0.0, a=:f)]
    @test_broken area(pm3d) ≈ 0.5 # Currently broken as zero(PointMeta(0.0, 0.0, 0.0, a=:d)) fails
end

@testset "embedding metadata" begin
    @testset "Meshes" begin

        @testset "per vertex attributes" begin
            points = rand(Point{3, Float64}, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(SVector{3, Float64}, 8)
            stress = LinRange(0, 1, 8)
            mesh = Mesh(meta(points, normals = normals, stress = stress), tfaces)

            @test hasproperty(coordinates(mesh), :stress)
            @test hasproperty(coordinates(mesh), :normals)
            @test coordinates(mesh).stress === stress
            @test coordinates(mesh).normals === normals
            @test coordinates(mesh).normals === normals
            @test GeometryBasics.faces(mesh) === tfaces
            @test propertynames(coordinates(mesh)) == (:position, :normals, :stress)

        end

        @testset "per face attributes" begin

            # Construct a cube out of Quads
            points = Point{3, Float64}[
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                (2.0, 2.0, 0.0), (0.0, 2.0, 0.0),
                (0.0, 0.0, 12.0), (2.0, 0.0, 12.0),
                (2.0, 2.0, 12.0), (0.0, 2.0, 12.0)
            ]

            facets = QuadFace{Cint}[
                1:4,
                5:8,
                [1,5,6,2],
                [2,6,7,3],
                [3, 7, 8, 4],
                [4, 8, 5, 1]
            ]

            markers = Cint[-1, -2, 0, 0, 0, 0]
            # attach some additional information to our faces!
            mesh = Mesh(points, meta(facets, markers = markers))
            @test hasproperty(GeometryBasics.faces(mesh), :markers)
            # test with === to assert we're not doing any copies
            @test GeometryBasics.faces(mesh).markers === markers
            @test coordinates(mesh) === points
            @test metafree(GeometryBasics.faces(mesh)) === facets

        end

    end

    @testset "polygon with metadata" begin
        polys = [Polygon(rand(Point{2, Float32}, 20)) for i in 1:10]
        pnames = [randstring(4) for i in 1:10]
        numbers = LinRange(0.0, 1.0, 10)
        bin = rand(Bool, 10)
        # create a polygon
        poly = PolygonMeta(polys[1], name = pnames[1], value = numbers[1], category = bin[1])
        # create a MultiPolygon with the right type & meta information!
        multipoly = MultiPolygonMeta(polys, name = pnames, value = numbers, category = bin)
        @test multipoly isa AbstractVector
        @test poly isa GeometryBasics.AbstractPolygon

        @test GeometryBasics.getcolumn(poly, :name) == pnames[1]
        @test GeometryBasics.MetaFree(PolygonMeta) == Polygon

        @test GeometryBasics.getcolumn(multipoly, :name) == pnames
        @test GeometryBasics.MetaFree(MultiPolygonMeta) == MultiPolygon

        meta_p = meta(polys[1], boundingbox=Rect(0, 0, 2, 2))
        @test meta_p.boundingbox === Rect(0, 0, 2, 2)
        @test metafree(meta_p) === polys[1]
        attributes(meta_p) == Dict{Symbol, Any}(:boundingbox => meta_p.boundingbox,
                                                :polygon => polys[1])
    end
    @testset "point with metadata" begin
        p = Point(1.1, 2.2)
        @test p isa AbstractVector{Float64}
        pm = PointMeta(1.1, 2.2; a=1, b=2)
        p1 = Point(2.2, 3.6)
        p2 = [p, p1]
        @test coordinates(p2) == p2
        @test meta(pm) === (a=1, b=2)
        @test metafree(pm) === p
        @test propertynames(pm) == (:position, :a, :b)
        @test GeometryBasics.MetaFree(typeof(pm)) == Point{2,Float64}
        @test_broken zero(pm) == [0, 0]
    end

    @testset "MultiPoint with metadata" begin
        p = collect(Point{2, Float64}(x, x+1) for x in 1:5)
        @test p isa AbstractVector
        mpm = MultiPointMeta(p, a=1, b=2)
        @test coordinates(mpm) == mpm
        @test meta(mpm) === (a=1, b=2)
        @test metafree(mpm) == p
        @test propertynames(mpm) == (:points, :a, :b)
    end

    @testset "LineString with metadata" begin
        linestring = LineStringMeta(Point{2, Int}[(10, 10), (20, 20), (10, 40)], a = 1, b = 2)
        @test linestring isa AbstractVector
        @test meta(linestring) === (a = 1, b = 2)
        @test metafree(linestring) == linestring
        @test propertynames(linestring) == (:lines, :a, :b)
    end

    @testset "MultiLineString with metadata" begin
        linestring1 = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
        linestring2 = LineString(Point{2, Int}[(40, 40), (30, 30), (40, 20), (30, 10)])
        multilinestring = MultiLineString([linestring1, linestring2])
        multilinestringmeta = MultiLineStringMeta([linestring1, linestring2]; boundingbox = Rect(1.0, 1.0, 2.0, 2.0))
        @test multilinestringmeta isa AbstractVector
        @test meta(multilinestringmeta) === (boundingbox = Rect(1.0, 1.0, 2.0, 2.0),)
        @test metafree(multilinestringmeta) == multilinestring
        @test propertynames(multilinestringmeta) == (:linestrings, :boundingbox)
    end

    @testset "Mesh with metadata" begin
       m = triangle_mesh(Sphere(Point3f(0), 1))
       m_meta = MeshMeta(m; boundingbox=Rect(1.0, 1.0, 2.0, 2.0))
       @test meta(m_meta) === (boundingbox = Rect(1.0, 1.0, 2.0, 2.0),)
       @test metafree(m_meta) === m
       @test propertynames(m_meta) == (:mesh, :boundingbox)
   end
end

@testset "embedding MetaT" begin
    @testset "MetaT{Polygon}" begin
        polys = [Polygon(rand(Point{2, Float32}, 20)) for i in 1:10]
        multipol = MultiPolygon(polys)
        pnames = [randstring(4) for i in 1:10]
        numbers = LinRange(0.0, 1.0, 10)
        bin = rand(Bool, 10)
        # create a polygon
        poly = MetaT(polys[1], name = pnames[1], value = numbers[1], category = bin[1])
        # create a MultiPolygon with the right type & meta information!
        multipoly = MetaT(multipol, name = pnames, value = numbers, category = bin)
        @test multipoly isa MetaT
        @test poly isa MetaT

        @test GeometryBasics.getcolumn(poly, :name) == pnames[1]
        @test GeometryBasics.getcolumn(multipoly, :name) == pnames

        meta_p = MetaT(polys[1], boundingbox=Rect(0, 0, 2, 2))
        @test meta_p.boundingbox === Rect(0, 0, 2, 2)
        @test GeometryBasics.metafree(meta_p) == polys[1]
        @test GeometryBasics.metafree(poly) == polys[1]
        @test GeometryBasics.metafree(multipoly) == multipol
        @test GeometryBasics.meta(meta_p) == (boundingbox = GeometryBasics.HyperRectangle{2,Int64}([0, 0], [2, 2]),)
        @test GeometryBasics.meta(poly) == (name = pnames[1], value = 0.0, category = bin[1])
        @test GeometryBasics.meta(multipoly) == (name = pnames, value = numbers, category = bin)
    end

    #=
    So mesh works differently for MetaT
    since `MetaT{Point}` not subtyped to `AbstractPoint`
    =#

   @testset "MetaT{Mesh}" begin
        @testset "per vertex attributes" begin
            points = rand(Point{3, Float64}, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(SVector{3, Float64}, 8)
            stress = LinRange(0, 1, 8)
            mesh_nometa = Mesh(points, tfaces)
            mesh = MetaT(mesh_nometa, normals = normals, stress = stress)

            @test hasproperty(mesh, :stress)
            @test hasproperty(mesh, :normals)
            @test mesh.stress == stress
            @test mesh.normals == normals
            @test GeometryBasics.faces(mesh.main) == tfaces
            @test propertynames(mesh) == (:main, :normals, :stress)
        end
    end
end

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

@testset "constructors" begin
    @testset "Mesh" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})

        mesh = Mesh(points, [1,2,3])
        @test mesh == [Triangle(points...)]

        x = Point{3}(1.0)
        mesh = Mesh([x], [TriangleFace(1, 1, 1)])
        @test mesh == [Triangle(x, x, x)]

        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        faces = connect([1, 2, 3, 4], SimplexFace{4})
        mesh = Mesh(points, faces)
        @test mesh == [Tetrahedron(points...)]

        points = rand(Point3f, 8)
        tfaces = [GLTriangleFace(1, 2, 3), GLTriangleFace(5, 6, 7)]
        normals = rand(Vec3f, 8)
        uv = rand(Vec2f, 8)
        mesh = Mesh(points, tfaces)
        meshuv = Mesh(meta(points; uv=uv), tfaces)
        meshuvnormal = Mesh(meta(points; normals=normals, uv=uv), tfaces)

        @test mesh isa GLPlainMesh
        @test meshuv isa GLUVMesh3D
        @test meshuvnormal isa GLNormalUVMesh3D

        t = Tesselation(Rect2f(0, 0, 2, 2), (30, 30))
        m = GeometryBasics.mesh(t, pointtype=Point3f, facetype=QuadFace{Int})
        m2 = GeometryBasics.mesh(m, facetype=QuadFace{GLIndex})
        @test GeometryBasics.faces(m2) isa Vector{QuadFace{GLIndex}}
        @test GeometryBasics.coordinates(m2) isa Vector{Point3f}

    end
end

@testset "decompose/triangulation" begin
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
    @test normals(m_normal) isa Vector{Vec3f}
    primitive = Rect2(0, 0, 1, 1)
    m_normal = normal_mesh(primitive)
    @test normals(m_normal) isa Vector{Vec3f}
    primitive = Rect3(0, 0, 0, 1, 1, 1)
    m_normal = normal_mesh(primitive)
    @test normals(m_normal) isa Vector{Vec3f}

    points = decompose(Point2f, Circle(Point2f(0), 1))
    tmesh = triangle_mesh(points)
    @test normals(tmesh) == nothing

    m = GeometryBasics.mesh(Sphere(Point3f(0), 1))
    @test normals(m) == nothing
    m_normals = pointmeta(m, Normal())
    @test normals(m_normals) isa Vector{Vec3f}

    @test texturecoordinates(m) == nothing
    r2 = Rect2(0.0, 0.0, 1.0, 1.0)
    @test collect(texturecoordinates(r2)) == [(0.0, 1.0), (1.0, 1.0), (0.0, 0.0), (1.0, 0.0)]
    r3 = Rect3(0.0, 0.0, 1.0, 1.0, 2.0, 2.0)
    @test first(texturecoordinates(r3)) == Vec3(0, 0, 0)
    uv = decompose_uv(m)
    @test Rect(Point.(uv)) == Rect(0, 0, 1, 1)

    points = decompose(Point2f, Circle(Point2f(0), 1))
    m = GeometryBasics.mesh(points)
    @test coordinates(m) === points

    linestring = LineString(Point{2, Int}[(10, 10), (20, 20), (10, 40)])
    pts = Point{2, Int}[(10, 10), (20, 20), (10, 40)]
    linestring = LineString(pts)
    pts_decomp = decompose(Point{2, Int}, linestring)
    @test pts === pts_decomp

    pts_ext = Point{2, Int}[(5, 1), (3, 3), (4, 8), (1, 2), (5, 1)]
    ls_ext = LineString(pts_ext)
    pts_int1 = Point{2, Int}[(2, 2), (3, 8),(5, 6), (3, 4), (2, 2)]
    ls_int1 = LineString(pts_int1)
    pts_int2 = Point{2, Int}[(3, 2), (4, 5),(6, 1), (1, 4), (3, 2)]
    ls_int2 =  LineString(pts_int2)
    poly_ext = Polygon(ls_ext)
    poly_ext_int = Polygon(ls_ext, [ls_int1, ls_int2])
    @test decompose(Point{2, Int}, poly_ext) == pts_ext
    @test decompose(Point{2, Int}, poly_ext_int) == [pts_ext..., pts_int1..., pts_int2...]
end

@testset "mesh" begin
    primitive = Triangle(Point2f(0), Point2f(1), Point2f(1,0))
    m = GeometryBasics.mesh(primitive)
    @test length(faces(m)) == 1
end

@testset "convert mesh + meta" begin
    m = uv_normal_mesh(Circle(Point2f(0), 1f0))
    # for 2D primitives we dont actually calculate normals
    @test !hasproperty(m, :normals)
end

@testset "convert mesh + meta" begin
    m = uv_normal_mesh(Rect3f(Vec3f(-1), Vec3f(1, 2, 3)))
    m_normal = normal_mesh(m)
    # make sure we don't loose the uv
    @test hasproperty(m_normal, :uv)
    @test m == m_normal
    # Make sure we don't create any copies
    @test m.position === m_normal.position
    @test m.normals === m_normal.normals
    @test m.uv === m_normal.uv

    m = GeometryBasics.mesh(Rect3f(Vec3f(-1), Vec3f(1, 2, 3));
                            uv=Vec2{Float64}, normaltype=Vec3{Float64}, pointtype=Point3{Float64})
    m_normal = normal_mesh(m)
    @test hasproperty(m_normal, :uv)
    @test m.position !== m_normal.position
    @test m.normals !== m_normal.normals
    # uv stays untouched, since we don't specify the element type in normalmesh
    @test m.uv === m_normal.uv
end

@testset "modifying meta" begin
    xx = rand(10)
    points = rand(Point3f, 10)
    m = GeometryBasics.Mesh(meta(points, xx=xx), GLTriangleFace[(1,2,3), (3,4,5)])
    color = rand(10)
    m = pointmeta(m; color=color)

    @test hasproperty(m, :xx)
    @test hasproperty(m, :color)
    @test_throws ErrorException GeometryBasics.MetaType(Simplex)
    @test_throws ErrorException GeometryBasics.MetaFree(Simplex)


    @test m.xx === xx
    @test m.color === color

    m, colpopt = GeometryBasics.pop_pointmeta(m, :color)
    m, xxpopt = GeometryBasics.pop_pointmeta(m, :xx)

    @test propertynames(m) == (:position,)
    @test colpopt === color
    @test xxpopt === xx

    @testset "creating meta" begin
        x = Point3f[(1,3,4)]
        # no meta gets added, so should stay the same
        @test meta(x) === x
        @test meta(x, value=[1]).position === x
    end
    pos = Point2f[(10, 2)]
    m = Mesh(meta(pos, uv=[Vec2f(1, 1)]), [GLTriangleFace(1, 1, 1)])
    @test m.position === pos
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
    @test tmesh isa GLPlainMesh
    @test coordinates(tmesh) === decompose(Point3f, tmesh)

    nmesh = normal_mesh(m)
    @test nmesh isa GLNormalMesh
    @test metafree(coordinates(nmesh)) === decompose(Point3f, nmesh)
    @test normals(nmesh) === decompose_normals(nmesh)

    m = GeometryBasics.mesh(s, pointtype=Point3f)
    @test m isa Mesh{3, Float32}
    @test coordinates(m) isa Vector{Point3f}
    @test GeometryBasics.faces(m) isa Vector{GLTriangleFace}
end

@testset "MetaT and heterogeneous data" begin
    ls = [LineString([Point(i, (i+1)^2/6), Point(i*0.86,i+5), Point(i/3, i/7)]) for i in 1:10]
    mls = MultiLineString([LineString([Point(i+1, (i)^2/6), Point(i*0.75,i+8), Point(i/2.5, i/6.79)]) for i in 5:10])
    poly = Polygon(Point{2, Int}[(40, 40), (20, 45), (45, 30), (40, 40)])
    geom = [ls..., mls, poly]
    prop = Any[(country_states = "India$(i)", rainfall = (i*9)/2) for i in 1:11]
    push!(prop, (country_states = 12, rainfall = 1000))   # a pinch of heterogeneity

    feat = [MetaT(i, j) for (i,j) = zip(geom, prop)]
    sa = meta_table(feat)

    @test nameof(eltype(feat)) == :MetaT
    @test eltype(sa) === MetaT{Any,(:country_states, :rainfall),Tuple{Any,Float64}}
    @test propertynames(sa) === (:main, :country_states, :rainfall)
    @test getproperty(sa, :country_states) isa Array{Any}
    @test getproperty(sa, :main) == geom

    maintype, metanames, metatype = GeometryBasics.getnamestypes(typeof(feat[1]))
    @test (metanames, metatype) == ((:country_states, :rainfall), Tuple{String,Float64})


    @test StructArrays.createinstance(typeof(feat[1]), LineString([Point(1, (2)^2/6), Point(1*0.86,6), Point(1/3, 1/7)]), "Mumbai", 100) isa typeof(feat[1])

    @test Base.getindex(feat[1], 1) isa Line
    @test Base.size(feat[1]) == (2,)
end

@testset "StructArrays integration" begin
    pt = meta(Point(0.0, 0.0), color="red", alpha=0.1)
    @test StructArrays.component(pt, :position) == Point(0.0, 0.0)
    @test StructArrays.component(pt, :color) == "red"
    @test StructArrays.component(pt, :alpha) == 0.1
    @test StructArrays.staticschema(typeof(pt)) ==
        NamedTuple{(:position, :color, :alpha), Tuple{Point2{Float64}, String, Float64}}
    @test StructArrays.createinstance(typeof(pt), Point(0.0, 0.0), "red", 0.1) == pt

    s = StructArray([pt, pt])
    @test StructArrays.components(s) == (
        position = [Point(0.0, 0.0), Point(0.0, 0.0)],
        color = ["red", "red"],
        alpha = [0.1, 0.1]
    )
    s[2] = meta(Point(0.1, 0.1), color="blue", alpha=0.3)
    @test StructArrays.components(s) == (
        position = [Point(0.0, 0.0), Point(0.1, 0.1)],
        color = ["red", "blue"],
        alpha = [0.1, 0.3]
    )
end

@testset "Tests from GeometryTypes" begin
    include("geometrytypes.jl")
end

@testset "GeoInterface" begin
    include("geointerface.jl")
end

using Aqua
# Aqua tests
# Intervals brings a bunch of ambiquities unfortunately
Aqua.test_all(GeometryBasics; ambiguities=false)

end  # testset "GeometryBasics"
