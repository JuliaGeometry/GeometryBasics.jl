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
    function ConicalFrustum{T}(baseCenter,baseRadius,topCenter,topRadius) where T

        all(isfinite.(baseCenter)) ? nothing : throw(ArgumentError("The bottom circle center point of a conical frustum needs to have finite coordinates and be a number."))

        all(isfinite.(topCenter)) ? nothing : throw(ArgumentError("The top circle center point of a conical frustum needs to have finite coordinates and be a number."))

        baseRadius > 0 ? nothing : throw(ArgumentError("The base radius of a conical frustum needs to be positive."))

        isfinite(baseRadius) ? nothing : throw(ArgumentError("The base radius of a conical frustum needs to be finite."))

        topRadius > 0 ? nothing : throw(ArgumentError("The top radius of a conical frustum needs to be positive."))

        isfinite(topRadius) ? nothing : throw(ArgumentError("The top radius of a conical frustum needs to be finite."))

        new{T}(baseCenter,baseRadius,topCenter,topRadius)

    end # function

end # struct

"""
An external constructor for converting inputs of different eltypes to a common type.
"""
function ConicalFrustum(baseCenter::Point3{T1},baseRadius::T2,topCenter::Point3{T3},topRadius::T4) where {T1,T2,T3,T4}

    T = promote_type(T1,T2,T3,T4)

    ConicalFrustum{T}(baseCenter,baseRadius,topCenter,topRadius)

end # function

"""
An external convenience constructor for creating instances from arrays instead of having to exlicitly create points.
"""
function ConicalFrustum(baseCenter::AbstractArray{T1},baseRadius::T2,topCenter::AbstractArray{T3},topRadius::T4) where {T1,T2,T3,T4}

    baseCenterPoint = Point3{T1}(baseCenter)

    topCenterPoint = Point3{T3}(topCenter)

    ConicalFrustum(baseCenterPoint,baseRadius,topCenterPoint,topRadius)

end # function

"""
An external convenience constructor for creating instances from NTuples instead of having to exlicitly create points.
"""
function ConicalFrustum(baseCenter::NTuple{3,T1},baseRadius::T2,topCenter::NTuple{3,T3},topRadius::T4) where {T1,T2,T3,T4}

    baseCenterPoint = Point3{T1}(baseCenter)

    topCenterPoint = Point3{T3}(topCenter)

    ConicalFrustum(baseCenterPoint,baseRadius,topCenterPoint,topRadius)

end # function

"""
An external convenience constructor for creating instances from Tuples instead of having to exlicitly create points.
"""
function ConicalFrustum(baseCenter::Tuple{T1,T2,T3},baseRadius::T4,topCenter::Tuple{T4,T5,T6},topRadius::T7) where {T1,T2,T3,T4,T5,T6,T7}

    T = promote_type(T1,T2,T3,T4,T5,T6,T7)

    baseCenterPoint = Point3{T}(baseCenter)

    topCenterPoint = Point3{T}(topCenter)

    ConicalFrustum(baseCenterPoint,baseRadius,topCenterPoint,topRadius)

end # function
