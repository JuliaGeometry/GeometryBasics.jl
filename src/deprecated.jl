using Base: @deprecate_binding

# Types ...f0 renamed to ...f
@deprecate_binding Vecf0 Vecf
@deprecate_binding Pointf0 Pointf
for i in 1:4
    for T in [:Point, :Vec]
        oldname = Symbol("$T$(i)f0")
        newname = Symbol("$T$(i)f")
        @eval begin
            @deprecate_binding $oldname $newname
        end
    end
    oldname = Symbol("Mat$(i)f0")
    newname = Symbol("Mat$(i)f")
    @eval begin
        @deprecate_binding $oldname $newname
    end
end

# Rect types
@deprecate_binding Rect2D  Rect2
@deprecate_binding Rect3D  Rect3
@deprecate_binding FRect   Rectf
@deprecate_binding FRect2D Rect2f
@deprecate_binding FRect3D Rect3f
@deprecate_binding IRect   Recti
@deprecate_binding IRect2D Rect2i
@deprecate_binding IRect3D Rect3i
@deprecate_binding TRect   RectT
