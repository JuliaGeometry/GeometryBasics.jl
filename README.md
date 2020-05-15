
[![Build Status](https://api.travis-ci.org/JuliaGeometry/GeometryBasics.jl.svg?branch=master)](https://travis-ci.com/SimonDanisch/GeometryBasics.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/SimonDanisch/GeometryBasics.jl?svg=true)](https://ci.appveyor.com/project/SimonDanisch/GeometryBasics-jl)
[![Codecov](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/SimonDanisch/GeometryBasics.jl)

**Documentation**:
[![Docs - stable](https://img.shields.io/badge/docs-stable-lightgrey.svg)](http://juliageometry.github.io/GeometryBasics.jl/stable/)
[![Docs - dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliageometry.github.io/GeometryBasics.jl/dev)

# GeometryBasics.jl

Basic Geometry Types.
This package aimes to offer a standard set of Geometry types, which easily work with metadata, query frameworks on geometries and different memory layouts.
The aim is to create a solid basis for Graphics/Plotting, finite elements analysis, Geo applications, and general geometry manipulations - while offering a julian API, that still allows performant C-interop.


It is planned to merge this package with [GeometryTypes](https://github.com/JuliaGeometry/GeometryTypes.jl/) once it's fully mature, or simply replace it!

To get an idea, look at the runtests:

```julia
using GeometryBasics
using GeometryBasics: Polygon, MultiPolygon, Point, LineFace, Polytope, Line
using GeometryBasics: Simplex, connect, Triangle, NSimplex, Tetrahedron
using GeometryBasics: QuadFace, hascolumn, getcolumn, metafree, coordinates, TetrahedronFace
using GeometryBasics: TupleView, TriangleFace, SimplexFace, LineString, Mesh, meta
using Test, Random, Query, StructArrays, Tables
using StaticArrays


@testset "embedding metadata" begin
    @testset "Meshes" begin

        @testset "per vertex attributes" begin
            points = rand(Point{3, Float64}, 8)
            tfaces = TetrahedronFace{Int}[(1, 2, 3, 4), (5, 6, 7, 8)]
            normals = rand(SVector{3, Float64}, 8)
            stress = LinRange(0, 1, 8)
            mesh = Mesh(meta(points, normals = normals, stress = stress), tfaces)

            @test hascolumn(coordinates(mesh), :stress)
            @test hascolumn(coordinates(mesh), :normals)
            @test coordinates(mesh).stress === stress
            @test coordinates(mesh).normals === normals
            @test coordinates(mesh).normals === normals
            @test GeometryBasics.faces(mesh) === tfaces

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
                [1, 5, 6, 2],
                [2, 6, 7, 3],
                [3, 7, 8, 4],
                [4, 8, 5, 1]
            ]

            markers = Cint[-1, -2, 0, 0, 0, 0]
            # attach some additional information to our faces!
            mesh = Mesh(points, meta(facets, markers = markers))
            @test hascolumn(GeometryBasics.faces(mesh), :markers)
            # test with === to assert we're not doing any copies
            @test getcolumn(GeometryBasics.faces(mesh), :markers) === markers
            @test coordinates(mesh) === points
            @test metafree(GeometryBasics.faces(mesh)) === facets

        end

    end
    @testset "polygon with metadata" begin
        polys = [Polygon(rand(Point{2, Float32}, 20)) for i in 1:10]
        pnames = [randstring(4) for i in 1:10]
        numbers = LinRange(0.0, 1.0, 10)
        bin = rand(Bool, 10)
        # create just an array
        plain = meta(polys, name = pnames, value = numbers, category = bin)
        # create a MultiPolygon with the right type & meta information!
        multipoly = MultiPolygon(polys; name = pnames, value = numbers, category = bin)
        for x in (plain, multipoly)
            for (mp, p, n, num, b) in zip(x, polys, pnames, numbers, bin)
                @test mp.polygon == p
                @test mp.name == n
                @test mp.value == num
                @test mp.category == b
            end

            filtered = @from i in x begin
                @where i.value < 0.7
                @select i
                @collect
            end
            @test length(filtered) == 7
        end
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

    @testset "face views" begin
        numbers = [1, 2, 3, 4, 5, 6]
        points = connect(numbers, Point{2})
        faces = connect([1, 2, 3], TriangleFace)
        triangles = connect(points, faces)
        @test triangles == [Triangle(Point(1, 2), Point(3, 4), Point(5, 6))]
        x = Point{3}(1.0)
        triangles = connect([x], [TriangleFace(1, 1, 1)])
        @test triangles == [Triangle(x, x, x)]
        points = connect([1, 2, 3, 4, 5, 6, 7, 8], Point{2})
        faces = connect([1, 2, 3, 4], SimplexFace{4})
        triangles = connect(points, faces)
        @test triangles == [Tetrahedron(points...)]
    end

end


@testset "constructors" begin
    @testset "LineFace" begin

        points = connect([1, 2, 3, 4, 5, 6], Point{2})
        linestring = LineString(points)
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        points = rand(Point{2, Float64}, 4)
        linestring = LineString(points, 2)
        @test linestring == [Line(points[1], points[2]), Line(points[3], points[4])]

        linestring = LineString([points[1] => points[2], points[2] => points[3]])
        @test linestring == [Line(points[1], points[2]), Line(points[2], points[3])]

        faces = [1, 2, 3]
        linestring = LineString(points, faces)
        @test linestring == LineString([points[1] => points[2], points[2] => points[3]])
        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        linestring = LineString(points, faces, 2)
        @test linestring == LineString([a => b, c => d])

        faces = [LineFace(1, 2), LineFace(3, 4)]
        linestring = LineString(points, faces)
        @test linestring == LineString([a => b, c => d])
    end

    @testset "Polygon" begin

        points = connect([1, 2, 3, 4, 5, 6], Point{2})
        polygon = Polygon(points)
        @test polygon == Polygon(LineString(points))

        points = rand(Point{2, Float64}, 4)
        linestring = LineString(points, 2)
        @test Polygon(points, 2) == Polygon(linestring)

        faces = [1, 2, 3]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))

        a, b, c, d = Point(1, 2), Point(3, 4), Point(5, 6), Point(7, 8)
        points = [a, b, c, d]; faces = [1, 2, 3, 4]
        polygon = Polygon(points, faces, 2)
        @test polygon == Polygon(LineString(points, faces, 2))

        faces = [LineFace(1, 2), LineFace(3, 4)]
        polygon = Polygon(points, faces)
        @test polygon == Polygon(LineString(points, faces))
    end

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

    end

end
```
