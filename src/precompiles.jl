using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        # Hits FaceView, QuadFace, all standard decompose's, some type conversions
        r = Rect3d(Point3d(0), Vec3d(1))
        m1 = uv_normal_mesh(r)

        c = Circle(Point2f(0), 1)
        m2 = uv_normal_mesh(c, pointtype = Point3f) # hits normal gen

        m = merge([m1, m2]) # hits mixed path, clear_faceviews, then normal path
        GeometryBasics.split_mesh(m)
        Rect3d(m)

        # Getters
        vertex_attributes(m)
        coordinates(m)
        normals(m)
        texturecoordinates(m)
        faces(m)
        
        face_normals(coordinates(r), faces(r))

        # Triangulation
        triangle_mesh(Polygon(rand(Point2f, 4)))

        # Other primitives
        uv_normal_mesh(Rect2(0,0,1,1))
        uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1f0), 3))
        uv_normal_mesh(Cylinder(Point3f(0), Point3f(0,0,1), 1f0))
        uv_normal_mesh(Pyramid(Point3f(0), 1f0, 1f0))

        # other disconnected compiles
        M = Mat3f(I)
        inv(M)
        M[1, Vec(1, 3)]
        M * Vec(1,2,3)

        Point2f(0.5, 0.1) in Triangle(Point2f(0), Point2f(0.5, 1), Point2f(1, 0))
        decompose(GLTriangleFace, [Point2f(0), Point2f(0.5, 1), Point2f(1, 0)])
        Point3f(0.5, 0, 1f0) in r 
    end
end
