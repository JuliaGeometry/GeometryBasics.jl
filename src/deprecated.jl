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
