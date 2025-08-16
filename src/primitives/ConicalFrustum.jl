"""
A type that represents a conical frustum: a cone whose top has been cut off such,
that the intersection of the cutting plane and the cone is a circle.
"""
struct ConicalFrustum{T} <: GeometryPrimitive{3,T}

    """
    The center point of the circular base.
    """
    baseCenter :: Point3{T}
    """
    The radius of the bottom circle.
    """
    baseRadius :: T

    """
    The center of the top circle.
    """
    topCenter :: Point3{T}

    """
    The radius of the top circle.
    """
    topRadius :: T

    """
    An inner constructor that validates the input values.
    """
    function ConicalFrustum(baseCenter,baseRadius,topCenter,topRadius)

        baseRadius > 0 ? nothing : throw(ArgumentError("The base radius of a conical frustum needs to be positive."))

        isfinite(baseRadius) ? nothing : throw(ArgumentError("The base radius of a conical frustum needs to be finite."))

        topRadius > 0 ? nothing : throw(ArgumentError("The top radius of a conical frustum needs to be positive."))

        isfinite(topRadius) ? nothing : throw(ArgumentError("The top radius of a conical frustum needs to be finite."))

        new(baseCenter,baseRadius,topCenter,topRadius)

    end # function

end # struct

"""
An external constructor for converting inputs of different eltypes to a common type.
"""
function ConicalFrustum(baseCenter::Point3{T1},baseRadius::T2,topCenter::Point3{T3},topRadius::T4) where {T1,T2,T3,T4}

    T = promote_type(T1,T2,T3,T4)

    ConicalFrustum{T}(baseCenter,baseRadius,topCenter,topRadius)

end # function
