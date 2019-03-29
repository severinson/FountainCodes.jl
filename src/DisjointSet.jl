############################################################
#
#   A forest of disjoint sets of integers

#   Extends the DataStructures.jl implementation by tracking the
#
#   number of vertices and edges in each component.
#
############################################################

export IntDisjointSetsTracked
export vertices, edges, cvertices, cedges, croot, reset!

mutable struct IntDisjointSetsTracked
    s::IntDisjointSets
    vertices::Vector{Int}
    edges::Vector{Int}
    maxv::Int
    maxe::Int
    maxroot::Int

    # creates a tracked disjoint set comprised of n singletons
    IntDisjointSetsTracked(n::Integer) = new(
        IntDisjointSets(n),
        ones(Int, n),
        zeros(Int, n),
        1, 0, 1,
    )
end

DataStructures.length(st::IntDisjointSetsTracked) = length(st.s)
DataStructures.num_groups(st::IntDisjointSetsTracked) = num_groups(st.s)
DataStructures.find_root(st::IntDisjointSetsTracked, x::Integer) = find_root(st.s, x)
DataStructures.in_same_set(st::IntDisjointSetsTracked, x::Integer, y::Integer) = in_same_set(st.s, x, y)

# find the number of edges/vertices in a component
vertices(st::IntDisjointSetsTracked, x::Integer) = st.vertices[find_root(st, x)]
edges(st::IntDisjointSetsTracked, x::Integer) = st.edges[find_root(st, x)]
cvertices(st::IntDisjointSetsTracked) = st.maxv
cedges(st::IntDisjointSetsTracked) = st.maxe
croot(st::IntDisjointSetsTracked) = st.maxroot

# merge the subset containing x and that containing y into one
# and return the root of the new set.
function DataStructures.union!(st::IntDisjointSetsTracked, x::Integer, y::Integer)
    v, e = vertices(st, x), edges(st, x)+1

    # TODO: can be made faster by writing a custom union! function
    # since union!(s, x, y) also compares vertices.
    if !in_same_set(st.s, x, y)
        v += vertices(st, y)
        e += edges(st, y)
    end
    root = union!(st.s, x, y)
    st.vertices[root], st.edges[root] = v, e
    if v >= st.maxv # store the root of the largest component
        st.maxv, st.maxe = v, e
        st.maxroot = root
    end
    root
end

# clear the data structure
function reset!(st::IntDisjointSetsTracked)
    @views st.s.parents[:] = 1:length(st.s.parents)
    @views st.s.ranks[:] = 0
    @views st.vertices[:] = 1
    @views st.edges[:] = 0
    st.maxv, st.maxe, st.maxroot = 1, 0, 1
    return
end
