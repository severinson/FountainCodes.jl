doc"Coded symbol representation."
struct R10Symbol
    esi::Int # encoded symbol id
    value::Int # value of the symbol
    neighbours::Set{Int} # indices of neighbouring symbols
    function R10Symbol(esi::Int, value::Int)
        new(esi, value, Set{Int}())
    end
    function R10Symbol(esi::Int, value::Int, neighbours::Set{Int})
        new(esi, value, neighbours)
    end
end
