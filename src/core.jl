# comaptibility with `GeometryBasicsCore` (dependencies)

Tables.schema(mesh::Mesh) = Tables.schema(getfield(mesh, :simplices))
Tables.schema(faceview::FaceView) = Tables.schema(getfield(faceview, :elements))

@inline function GeometryBasicsCore.connect(points::AbstractMatrix{T},
                         P::Type{<:AbstractPoint{N}}) where {T <: Real,N}
    return if size(points, 1) === N
        return reinterpret(Point{N,T}, points)
    elseif size(points, 2) === N
        seglen = size(points, 1)
        columns = ntuple(N) do i
            return view(points, ((i - 1) * seglen + 1):(i * seglen))
        end
        return StructArray{Point{N,T}}(columns)
    else
        error("Dim 1 or 2 must be equal to the point dimension!")
    end
end
