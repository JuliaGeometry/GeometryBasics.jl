function decompose_triangulate_fallback(primitive::Meshable, pointtype, facetype)
    positions = decompose(pointtype, primitive)
    faces = decompose(facetype, primitive)
    # If faces returns nothing for primitive, we try to triangulate!
    if faces === nothing
        # triangulation.jl
        faces = decompose(facetype, positions)
    end
    return positions, faces
end

"""
    mesh(primitive::GeometryPrimitive;
         pointtype=Point, facetype=GLTriangle,
         uv=nothing, normaltype=nothing)

Creates a mesh from `primitive`.

Uses the element types from the keyword arguments to create the attributes.
The attributes that have their type set to nothing are not added to the mesh.
Note, that this can be an `Int` or `Tuple{Int, Int}``, when the primitive is grid based.
It also only losely correlates to the number of vertices, depending on the algorithm used.
#TODO: find a better number here!
"""
function mesh(primitive::Meshable; pointtype=Point, facetype=GLTriangleFace)
    positions, faces = decompose_triangulate_fallback(primitive, pointtype, facetype)
    return Mesh(positions, faces)
end

"""
    mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
         normaltype=nothing)

Create a mesh from a polygon given as a vector of points, using triangulation.
"""
function mesh(polygon::AbstractVector{P}; pointtype=P, facetype=GLTriangleFace,
              normaltype=nothing) where {P<:Point{2}}

    return mesh(Polygon(polygon); pointtype=pointtype, facetype=facetype,
                normaltype=normaltype)
end

function mesh(polygon::AbstractPolygon{Dim,T}) where {Dim,T}
    positions, faces = decompose_triangulate_fallback(polygon; pointtype=Point{Dim, T}, facetype=GLTriangleFace)
    return Mesh(positions, faces)
end

function triangle_mesh(primitive::Meshable{N}) where {N}
    return Mesh(decompose(Point{N,Float32}, primitive), decompose(GLTriangleFace, primitive))
end

function uv_mesh(primitive::Meshable{N,T}) where {N,T}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (uv=decompose_uv(primitive),))
end

function uv_normal_mesh(primitive::Meshable{N}) where {N}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (uv=decompose_uv(primitive), normals=decompose_normals(primitive)))
end

function normal_mesh(points::AbstractVector{<:Point},
                     faces::AbstractVector{<:AbstractFace})
    _points = decompose(Point3f, points)
    _faces = decompose(GLTriangleFace, faces)
    return MetaMesh(Mesh(_points, _faces), (normals=normals(_points, _faces),))
end

function normal_mesh(primitive::Meshable{N}) where {N}
    m = triangle_mesh(primitive)
    return MetaMesh(m, (normals=decompose_normals(primitive),))
end

"""
    volume(triangle)

Calculate the signed volume of one tetrahedron. Be sure the orientation of your
surface is right.
"""
function volume(triangle::Triangle) where {VT,FT}
    v1, v2, v3 = triangle
    sig = sign(orthogonal_vector(v1, v2, v3) ⋅ v1)
    return sig * abs(v1 ⋅ (v2 × v3)) / 6
end

"""
    volume(mesh)

Calculate the signed volume of all tetrahedra. Be sure the orientation of your
surface is right.
"""
function volume(mesh::Mesh) where {VT,FT}
    return sum(volume, mesh)
end

function Base.merge(meshes::AbstractVector{<:Mesh})
    return if isempty(meshes)
        error("No meshes to merge")
    elseif length(meshes) == 1
        return meshes[1]
    else
        m1 = meshes[1]
        ps = copy(coordinates(m1))
        fs = copy(faces(m1))
        for mesh in Iterators.drop(meshes, 1)
            append!(fs, map(f -> f .+ length(ps), faces(mesh)))
            append!(ps, coordinates(mesh))
        end
        return Mesh(ps, fs)
    end
end
