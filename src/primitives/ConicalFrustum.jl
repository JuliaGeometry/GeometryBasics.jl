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

# Accessor functions for frustum fields.

"""
An accessor function for a base center.
"""
baseCenter(x::ConicalFrustum) = x.baseCenter

"""
An accessor function for base radius.
"""
baseRadius(x::ConicalFrustum) = x.baseRadius

"""
An accessor function for top center.
"""
topCenter(x::ConicalFrustum) = x.topCenter

"""
An accessor function for top radius.
"""
topRadius(x::ConicalFrustum) = x.topRadius

# Functions for computing derived properties based on type field values.

"""
Computes the length of a frustum as the norm of the difference of its top and base centers.
"""
Base.length(x::ConicalFrustum) = LinearAlgebra.norm(topCenter(x) - baseCenter(x))

"""
Computes the area of the base of a given conical frustum.
"""
baseArea(x::ConicalFrustum) = pi * baseRadius(x) ^ 2

"""
Computes the area of the top of a given conical frustum.
"""
topArea(x::ConicalFrustum) = pi * topRadius(x) ^ 2

"""
Computes the slant length of a conical frustum.
"""
slantLength(x::ConicalFrustum) = sqrt( ( baseRadius(x) - topRadius(x) ) ^ 2 + length(x) ^ 2 )

"""
Computes the surface area of a frustum not including the top and bottom areas.
"""
function slantArea(x::ConicalFrustum)

    baseRadiusVal = baseRadius(x)

    topRadiusVal = topRadius(x)

    lengthVal = length(x)

    pi * (baseRadiusVal + topRadiusVal) * slantLength(x)

end # function

"""
Computes the total surface area of a conical frustum.
"""
surfaceArea(x::ConicalFrustum) = baseArea(x) + topArea(x) + slantArea(x)

"""
Computes the volume of a conical frustum.
"""
function volume(x::ConicalFrustum)

    baseRadiusVal = baseRadius(x)

    topRadiusVal = topRadius(x)

    lengthVal = length(x)

    pi * lengthVal * ( baseRadiusVal ^ 2 + baseRadiusVal * topRadiusVal + topRadiusVal ^ 2 ) / 3

end # function

"""
Computes the centroid of a frustum: the mean between the base and top centers.
"""
centroid(x::ConicalFrustum) = (baseCenter(x) + topCenter(x)) / 2

"""
Computes the coordinates required for the discretization
of a frustum. The logic is the same as that for a cylinder,
where the top and bottom circles are approximated using
a polygon.
"""
function coordinates(c::ConicalFrustum{T}, nvertices=30) where {T}

    nvertices += isodd(nvertices)

    nhalf = div(nvertices, 2)

    R = rotation(c)

    step = 2pi / nhalf

    ps = Vector{Point3{T}}(undef, nvertices + 2)

    baseRadiusVal = baseRadius(c)

    baseCenterVal = baseCenter(c)

    topRadiusVal = topRadius(c)

    topCenterVal = topCenter(c)

    # First discretize the base...

    for i in 1:nhalf

        phi = (i-1) * step

        ps[i] = R * Point3{T}(baseRadius * cos(phi), baseRadius * sin(phi), 0) + baseCenterVal

    end

    # ... and then the top circle.

    for i in 1:nhalf

        phi = (i-1) * step

        ps[i + nhalf] = R * Point3{T}(topRadiusVal * cos(phi), topRadiusVal * sin(phi), 0) + topCenterVal
    end

    ps[end-1] = baseCenterVal

    ps[end] = topCenterVal

    return ps

end # function
