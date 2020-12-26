function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    @assert precompile(HyperRectangle{2,Float32}, (Int, Int, Int, Int))
    @assert precompile(==, (HyperRectangle{2,Float32}, HyperRectangle{2,Float32}))
    @assert precompile(normal_mesh, (Tesselation{3,Float32,Cylinder{3,Float32},1},))
    @assert precompile(normal_mesh, (Tesselation{3,Float32,HyperSphere{3,Float32},1},))
    @assert precompile(normal_mesh, (HyperSphere{3,Float32},))
    @assert precompile(Tuple{typeof(*),SMatrix{4, 4, Float32, 16},HyperRectangle{3, Float32}})   # time: 0.11091917
    @assert precompile(Tuple{typeof(coordinates),HyperRectangle{2, Float32},Tuple{Int64, Int64}})   # time: 0.08693867
    @assert precompile(union, (HyperRectangle{3, Float32}, HyperRectangle{3, Float32}))
    @assert precompile(Tuple{typeof(decompose),Type{Point{2, Float32}},HyperRectangle{2, Float32}})   # time: 0.026609203
    @assert precompile(Tuple{Type{HyperRectangle{3, Float32}},HyperRectangle{2, Float32}})   # time: 0.023717888
    @assert precompile(Tuple{typeof(+),HyperRectangle{3, Float32},Point{3, Float32}})   # time: 0.006633118
    @assert precompile(Tuple{Type{Rect2D{T} where T},Float32,Float32,Float32,Float32})   # time: 0.001636267
    @assert precompile(Tuple{typeof(*),HyperRectangle{2, Float32},Float32})   # time: 0.001057589

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
