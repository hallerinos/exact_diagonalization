function buildAdjacencyMatrix(L::Array{Int64,1})
    size = prod(L)
    # build the adjacency matrix
    AM = zeros(size, size)
    for r=1:size
        mod(r,L[1])!=0 && r+1<=size ? AM[r, r+1] = 1 : Nothing
        r+L[1]<=size ? AM[r, r+L[1]] = 1 : Nothing
        mod(r,L[1])!=0 && r+L[1]+1<=size ? AM[r, r+L[1]+1] = 1 : Nothing
    end
    # show(spy(AM))
    # display(gplot(DiGraph(AM), nodelabel=1:size))
    return AM
end