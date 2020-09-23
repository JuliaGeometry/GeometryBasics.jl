#! /usr/bin/env julia

using Documenter

using GeometryBasics

# Copy the README to serve as the homepage
cp(joinpath(@__DIR__, "..", "README.md"), joinpath(@__DIR__, "src", "index.md"))

makedocs(format=Documenter.HTML(), sitename="GeometryBasics.jl",
         pages=["index.md", "primitives.md", "rectangles.md", "polygons.md", "meshes.md",
                "decomposition.md", "distancefields.md", "metadata.md", "api.md"],
         modules=[GeometryBasics])

deploydocs(repo="github.com/JuliaGeometry/GeometryBasics.jl.git", push_preview=true)
