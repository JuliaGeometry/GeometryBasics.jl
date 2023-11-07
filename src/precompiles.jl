
function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Point2f(0.5, 0.1) in Triangle(Point2f(0), Point2f(0.5, 1), Point2f(1, 0))
    decompose(GLTriangleFace, [Point2f(0), Point2f(0.5, 1), Point2f(1, 0)])
end
