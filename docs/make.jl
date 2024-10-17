#! /usr/bin/env julia

using Documenter
using GeometryBasics

DocMeta.setdocmeta!(GeometryBasics, :DocTestSetup, :(using GeometryBasics); recursive=true)

makedocs(format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true"),
         sitename="GeometryBasics.jl",
         pages=[
                "index.md",
                "primitives.md",
                "polygons.md",
                "meshes.md",
                "decomposition.md",
                "static_array_types.md",
                "api.md"
               ],
         modules=[GeometryBasics])

deploydocs(repo="github.com/JuliaGeometry/GeometryBasics.jl.git", push_preview=true)
