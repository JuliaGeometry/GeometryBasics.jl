module GeometryBasics

    using StaticArrays, Tables, StructArrays, IterTools, LinearAlgebra

    using Base: @propagate_inbounds

    include("fixed_arrays.jl")
    include("offsetintegers.jl")
    include("basic_types.jl")
    include("metadata.jl")
    include("viewtypes.jl")
    include("geometry_primitives.jl")
    include("rectangles.jl")
    include("triangulation.jl")
    include("meshes.jl")

end # module
