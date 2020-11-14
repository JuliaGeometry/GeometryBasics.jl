function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    @assert precompile(HyperRectangle{2,Float32}, (Int, Int, Int, Int))
    @assert precompile(==, (HyperRectangle{2,Float32}, HyperRectangle{2,Float32}))
    @assert precompile(normal_mesh, (Tesselation{3,Float32,Cylinder{3,Float32},1},))
    @assert precompile(normal_mesh, (Tesselation{3,Float32,HyperSphere{3,Float32},1},))
    @assert precompile(normal_mesh, (HyperSphere{3,Float32},))

    if Base.VERSION >= v"1.6.0-DEV.1083"
        @assert precompile(triangle_mesh, (Polygon{2, Float32, Point2f0, LineString{2, Float32, Point2f0,
                                        Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f0, Point2f0}, TupleView{Tuple{Point2f0, Point2f0}, 2, 1, Vector{Point2f0}}, false}},
                                        Vector{LineString{2, Float32, Point2f0, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f0, Point2f0}, TupleView{Tuple{Point2f0, Point2f0}, 2, 1, Vector{Point2f0}}, false}}}},))
    else
        @assert precompile(triangle_mesh, (Polygon{2, Float32, Point2f0, LineString{2, Float32, Point2f0,
                                        Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f0, Point2f0}, TupleView{Tuple{Point2f0, Point2f0}, 2, 1, Vector{Point2f0}}}},
                                        Vector{LineString{2, Float32, Point2f0, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f0, Point2f0}, TupleView{Tuple{Point2f0, Point2f0}, 2, 1, Vector{Point2f0}}}}}},))
    end

    @assert precompile(split_intersections, (Vector{Point2f0},))
end
