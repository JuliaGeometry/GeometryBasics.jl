module GeometryBasics

using StaticArrays, Tables, StructArrays
using Requires

using Base: @propagate_inbounds

include("FixedSizeArrays.jl")
using .FixedSizeArrays

include("basic_types.jl")
include("metadata.jl")
include("viewtypes.jl")


function __init__()
    @require AbstractPlotting = "537997a7-5e4e-5d89-9595-2241ea00577e" include("abstractplotting.jl")
end

end # module
