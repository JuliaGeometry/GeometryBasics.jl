#! /usr/bin/env julia

using Documenter
using GeometryBasics

DocMeta.setdocmeta!(GeometryBasics, :DocTestSetup, :(using GeometryBasics); recursive=true)

makedocs(; sitename="GeometryBasics.jl",
         format=Documenter.HTML(; prettyurls=false, size_threshold=3000000,
                                example_size_threshold=3000000),
         pages=["index.md",
                "primitives.md",
                "polygons.md",
                "meshes.md",
                "decomposition.md",
                "boundingboxes.md",
                "static_array_types.md",
                "api.md"],
         modules=[GeometryBasics])

deploydocs(; repo="github.com/JuliaGeometry/GeometryBasics.jl.git", push_preview=true)
