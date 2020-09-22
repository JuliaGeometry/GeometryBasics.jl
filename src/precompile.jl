function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    precompile(normal_mesh, (Tesselation{3,Float32,Cylinder{3,Float32},1},))
    precompile(normal_mesh, (Tesselation{3,Float32,HyperSphere{3,Float32},1},))
    precompile(normal_mesh, (HyperSphere{3,Float32},))
end
