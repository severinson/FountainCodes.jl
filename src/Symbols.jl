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
    active_neighbours::Array{Int,1}
    inactive_neighbours::Array{Int,1}
    function R10Symbol(esi::Int, value::Int,
                       primary_neighbour::Int,
                       active_neighbours::Array{Int,1},
                       inactive_neighbours::Array{Int,1},
                       sort=true)
        if sort
            return new(
                esi,
                value,
                primary_neighbour,
                sort!(copy(active_neighbours)),
                sort!(copy(inactive_neighbours)),
            )
        else
            return new(
                esi,
                value,
                primary_neighbour,
                active_neighbours,
                inactive_neighbours,
            )
        end
    end
    function R10Symbol(esi::Int, value::Int, neighbours::Array{Int,1})
        R10Symbol(esi, value, -1, neighbours, Array{Int,1}())
    end
end

doc"Number of neighbouring outer coded symbols."
function degree(is::ISymbol)
    return length(is.neighbours)
end

doc"Number of neighbouring intermediate symbols."
function degree(cs::R10Symbol)
    return length(cs.active_neighbours) + length(cs.inactive_neighbours)
end

doc"Neighbouring outer code symbols."
function neighbours(is::ISymbol)
    return collect(is.neighbours)
end

doc"Neighbouring intermediate symbols."
function neighbours(cs::R10Symbol)
    return append!(append!([], cs.active_neighbours), cs.inactive_neighbours)
end

doc"Number of non-zero entries in V."
function active_degree(cs::R10Symbol)
    return length(cs.active_neighbours)
end

doc"Number of non-zero entries in U."
function inactive_degree(cs::R10Symbol)
    return length(cs.inactive_neighbours)
end

doc"Neighbours that are not decoded or inactivated."
function active_neighbours(cs::R10Symbol)
    return cs.active_neighbours
end
