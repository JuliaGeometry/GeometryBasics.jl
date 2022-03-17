m = GeometryBasics.normal_mesh(Sphere(Point3f(0), 1f0))
@test decompose_uv(m) isa Vector{Vec2f}

o, extr, r = Point2f(1, 2), Point2f(3, 4), 5.0f0
s = Cylinder(o, extr, r)
positions = Point{3,Float32}[(-0.7677671, 3.767767, 0.0),
                                (2.767767, 0.23223293, 0.0),
                                (0.23223293, 4.767767, 0.0),
                                (3.767767, 1.2322329, 0.0), (1.2322329, 5.767767, 0.0),
                                (4.767767, 2.232233, 0.0)]
@test decompose(Point3f, Tesselation(s, (2, 3))) ≈ positions

FT = TriangleFace{Int}
f = FT[(1, 2, 4), (1, 4, 3), (3, 4, 6), (3, 6, 5)]
@test f == decompose(FT, Tesselation(s, (2, 3)))

v1 = Point{3,Float64}(1, 2, 3)
v2 = Point{3,Float64}(4, 5, 6)
R = 5.0
s = Cylinder(v1, v2, R)
positions = Point{3,Float64}[(4.535533905932738, -1.5355339059327373, 3.0),
                                (7.535533905932738, 1.4644660940672627, 6.0),
                                (3.0412414523193148, 4.041241452319315,
                                -1.0824829046386295),
                                (6.041241452319315, 7.041241452319315,
                                1.9175170953613705),
                                (-2.535533905932737, 5.535533905932738,
                                2.9999999999999996),
                                (0.46446609406726314, 8.535533905932738, 6.0),
                                (-1.0412414523193152, -0.04124145231931431,
                                7.0824829046386295),
                                (1.9587585476806848, 2.9587585476806857,
                                10.08248290463863), (1, 2, 3), (4, 5, 6)]

@test decompose(Point3{Float64}, Tesselation(s, 8)) ≈ positions

f = TriangleFace{Int}[(3, 2, 1), (4, 2, 3), (5, 4, 3), (6, 4, 5), (7, 6, 5),
                            (8, 6, 7), (1, 8, 7), (2, 8, 1), (3, 1, 9), (2, 4, 10),
                            (5, 3, 9), (4, 6, 10), (7, 5, 9), (6, 8, 10), (1, 7, 9),
                            (8, 2, 10)]
@test f == decompose(TriangleFace{Int}, Tesselation(s, 8))

m = triangle_mesh(Tesselation(s, 8))

@test GeometryBasics.faces(m) == f
