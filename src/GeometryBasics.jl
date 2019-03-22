module GeometryBasics

using StaticArrays, Tables, StructArrays

using Base: @propagate_inbounds

include("static_arrays.jl")
include("basic_types.jl")
include("metadata.jl")
include("viewtypes.jl")
include("simple.jl")
include("offsetinteger.jl")
include("hyperrectangles.jl")


end # module
