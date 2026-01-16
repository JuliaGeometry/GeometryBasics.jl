module ExtentsExt

import GeometryBasics
import Extents

function Extents.extent(rect::GeometryBasics.Rect2)
    (xmin, ymin), (xmax, ymax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end

function Extents.extent(rect::GeometryBasics.Rect3)
    (xmin, ymin, zmin), (xmax, ymax, zmax) = extrema(rect)
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax), Z=(zmin, zmax))
end

end
