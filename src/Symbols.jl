export BSymbol, QSymbol

"outer code symbol with binary coefficients."
struct BSymbol{VT} <: CodeSymbol
    esi::Int # encoded symbol id
    value::VT # value of the symbol
    neighbours::Vector{Int}
end

"number of neighbouring intermediate symbols."
function degree(cs::CodeSymbol)
    return length(cs.neighbours)
end

"outer code symbol with arbitrary coefficient type."
struct QSymbol{VT,CT} <: CodeSymbol
    esi::Int # encoded symbol id
    value::VT # value of the symbol
    neighbours::Vector{Int}
    coefficients::Vector{CT}
end
