
[![CI](https://github.com/JuliaGeometry/GeometryBasics.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaGeometry/GeometryBasics.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/GeometryBasics.jl)

[![Docs - stable](https://img.shields.io/badge/docs-stable-lightgrey.svg)](http://juliageometry.github.io/GeometryBasics.jl/stable/)
[![Docs - dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliageometry.github.io/GeometryBasics.jl/dev/)

# GeometryBasics.jl

Basic geometry types.

This package aims to offer a standard set of geometry types that easily work
with metadata, query frameworks on geometries and different memory layouts. The
aim is to create a solid basis for graphics/plotting, finite element analysis,
geo applications, and general geometry manipulations - while offering a Julian
API that still allows performant C-interop.

This package is a replacement for the discontinued [GeometryTypes](https://github.com/JuliaGeometry/GeometryTypes.jl/).

**Documentation:** http://juliageometry.github.io/GeometryBasics.jl/stable/

## Contributing

Make sure your changes don't break the documentation.

To build the documentation locally, you first need to instantiate the `docs/` project:

```
julia --project=docs/
pkg> instantiate
pkg> dev .
```

Then use `julia --project=docs/ docs/make.jl` to build the documentation. This
will also run the doctests defined in Markdown files. The doctests should be
written for the Julia version configured in [ci.yml](.github/workflows/ci.yml)
(`:docs` section).
