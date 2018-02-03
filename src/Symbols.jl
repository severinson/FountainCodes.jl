doc"Arbitrary coded symbol."
abstract type CodeSymbol end

doc"Intermediate code symbol."
struct ISymbol <: CodeSymbol
    value::Int
    neighbours::Set{Int}
    function ISymbol(value::Int, neighbours::Set{Int})
        new(value, neighbours)
    end
    function ISymbol(value::Int)
        new(value, Set{Int}())
    end
end

doc"Outer code symbol."
struct R10Symbol <: CodeSymbol
    esi::Int # encoded symbol id
    value::Int # value of the symbol
    primary_neighbour::Int
    active_neighbours::Set{Int}
    inactive_neighbours::Set{Int}
    # function R10Symbol(esi::Int, value::Int)
    #     new(esi, value, Set{Int}())
    # end
    function R10Symbol(esi::Int, value::Int,
                       primary_neighbour::Int,
                       active_neighbours::Set{Int},
                       inactive_neighbours::Set{Int})
        new(
            esi,
            value,
            primary_neighbour,
            active_neighbours,
            inactive_neighbours,
        )
    end
    function R10Symbol(esi::Int, value::Int, neighbours::Set{Int})
        R10Symbol(esi, value, -1, neighbours, Set{Int}())
    end
end
