macro warnpcfail(ex::Expr)
    modl = __module__
    file = __source__.file === nothing ? "?" : String(__source__.file)
    line = __source__.line
    quote
        $(esc(ex)) || @warn """precompile directive
     $($(Expr(:quote, ex)))
 failed. Please report an issue in $($modl) (after checking for duplicates) or remove this directive.""" _file=$file _line=$line
    end
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    @warnpcfail precompile(HyperRectangle{2,Float32}, (Int, Int, Int, Int))
    @warnpcfail precompile(==, (HyperRectangle{2,Float32}, HyperRectangle{2,Float32}))
    @warnpcfail precompile(Tuple{typeof(*),Mat{4, 4, Float32, 16},HyperRectangle{3, Float32}})   # time: 0.11091917
    @warnpcfail precompile(Tuple{typeof(coordinates),HyperRectangle{2, Float32},Tuple{Int64, Int64}})   # time: 0.08693867
    @warnpcfail precompile(union, (HyperRectangle{3, Float32}, HyperRectangle{3, Float32}))
    @warnpcfail precompile(Tuple{typeof(decompose),Type{Point{2, Float32}},HyperRectangle{2, Float32}})   # time: 0.026609203
    @warnpcfail precompile(Tuple{Type{HyperRectangle{3, Float32}},HyperRectangle{2, Float32}})   # time: 0.023717888
    @warnpcfail precompile(Tuple{typeof(+),HyperRectangle{3, Float32},Point{3, Float32}})   # time: 0.006633118
    @warnpcfail precompile(Tuple{Type{Rect2{T} where T},Float32,Float32,Float32,Float32})   # time: 0.001636267
    @warnpcfail precompile(Tuple{typeof(*),HyperRectangle{2, Float32},Float32})   # time: 0.001057589
    @warnpcfail precompile(split_intersections, (Vector{Point2f},))
end
