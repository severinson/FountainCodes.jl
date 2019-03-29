"""
    SimpleSelect

Benchmark selector.

"""
struct SimpleSelect <: Selector
    pq::PriorityQueue{Int,Int,Base.Order.ForwardOrdering}
    function SimpleSelect()
        new(PriorityQueue{Int,Int}())
    end
end

function Base.push!(sel::SimpleSelect, d::Decoder, rpi::Int, row::Row)
    enqueue!(sel.pq, rpi, degree(row))
end

function Base.pop!(sel::SimpleSelect, d::Decoder)
    while length(sel.pq) > 0
        rpi = dequeue!(sel.pq)
        row = d.rows[rpi]
        vdeg = vdegree(d, row)
        if vdeg > 0
            return d.rowperminv[rpi]
        end
    end
    error("row rows of non-zero vdegree")
end
