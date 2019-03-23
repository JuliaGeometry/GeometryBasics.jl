module GeometryBasics

using StaticArrays, Tables, StructArrays

using Base: @propagate_inbounds

include("FixedSizeArrays.jl")
using .FixedSizeArrays

include("basic_types.jl")
include("metadata.jl")
include("viewtypes.jl")
include("abstractplotting.jl")


end # module
