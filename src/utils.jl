const typesuffix = Dict(Float32=>"f", Int=>"i")

macro alias(T, dims, types...)
    for i in eval(dims)
        name = Symbol("$T$i")
        @eval begin
            const $name = $T{$i}
            export $name
        end
    end
    for t in eval(types)
        name = Symbol("$T$(typesuffix[eval(t)])")
        @eval begin
            const $name{N} = $T{N,$t}
            export $name
        end
    end
    for i in eval(dims)
        for t in eval(types)
            name = Symbol("$T$(i)$(typesuffix[eval(t)])")
            @eval begin
                const $name = $T{$i,$t}
                export $name
            end
        end
    end
end