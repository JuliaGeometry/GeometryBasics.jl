import AbstractPlotting
import AbstractPlotting.GeometryTypes

function to_triangle(triangles, tetra)
    push!(triangles, GeometryTypes.GLTriangle(tetra[1], tetra[2], tetra[3]))
    push!(triangles, GeometryTypes.GLTriangle(tetra[2], tetra[3], tetra[4]))
    push!(triangles, GeometryTypes.GLTriangle(tetra[2], tetra[3], tetra[4]))
    push!(triangles, GeometryTypes.GLTriangle(tetra[3], tetra[1], tetra[4]))
end
function to_lines(lines, tetra)
    push!(lines, tetra[1], tetra[2])
    push!(lines, tetra[2], tetra[3])
    push!(lines, tetra[3], tetra[4])
    push!(lines, tetra[2], tetra[4])
    push!(lines, tetra[1], tetra[4])
end


function AbstractPlotting.convert_arguments(T::Type{<: AbstractPlotting.Lines}, mesh::AbstractVector{<: GeometryBasics.Tetrahedron})
    points = GeometryTypes.Point3f0[]
    for elem in mesh
        to_lines(points, elem)
    end
    return (points,)
end

function AbstractPlotting.convert_arguments(::Type{<: AbstractPlotting.Mesh}, mesh::GeometryBasics.Mesh{T, N, <: GeometryBasics.Tetrahedron}) where {T, N}
    points = GeometryTypes.Point3f0.(coordinates(mesh))
    triangles = GeometryTypes.GLTriangle[]
    for elem in faces(mesh)
        to_triangle(triangles, elem)
    end
    return (GeometryTypes.GLNormalMesh(points, triangles),)
end

function AbstractPlotting.convert_arguments(T::Type{<: AbstractPlotting.Wireframe}, mesh::GeometryBasics.Mesh)
    m = convert_arguments(AbstractPlotting.Mesh, mesh)[1]
    return convert_arguments(T, m)
end
function AbstractPlotting.convert_arguments(p::AbstractPlotting.PointBased, mesh::GeometryBasics.Mesh)
    return convert_arguments(p, GeometryTypes.Point.(coordinates(mesh)))
end

function AbstractPlotting.convert_arguments(::Type{<: AbstractPlotting.Mesh}, mesh::GeometryBasics.Mesh)
    (GeometryTypes.GLNormalMesh(
        GeometryTypes.Point3f0.(coordinates(mesh)),
        GeometryTypes.GLTriangle.(faces(mesh))
    ),)
end


function to_mesh(mesh::GeometryTypes.AbstractMesh)
    points = GeometryBasics.Point{3, Float64}.(GeometryTypes.vertices(mesh))
    facets = GeometryBasics.TriangleFace{Cint}.(GeometryTypes.faces(mesh))
    Mesh(points, facets)
end
