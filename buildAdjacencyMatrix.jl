function buildAdjacencyMatrixRegular(L::Array{Int64,1})
    numsites = prod(L)
    # build the adjacency matrix
    AM = zeros(numsites, numsites)
    for r=1:numsites
        # the red links
        mod(r,L[1])!=0 && r+1<=numsites ? AM[r, r+1] = 1 : Nothing
        # the blue links
        mod(r,L[1])!=0 && r+L[1]+1<=numsites ? AM[r, r+L[1]+1] = 2 : Nothing
        # the green links
        r+L[1]<=numsites ? AM[r, r+L[1]] = 3 : Nothing
    end
    return AM,numsites
end

function buildAdjacencyMatrixSnowflake(L::Int64)
    numsites = 1 - 3*L + 3*L^2

    # build the adjacency matrix
    AM = zeros(numsites, numsites)
    # the red links
    for y=1:L-1
        AM[y,y+1] = 1
    end
    for y=L+1:2*L+1-1
        AM[y,y+1] = 1
    end
    for y=2*L+1+1:3*L+3-1
        AM[y,y+1] = 1
    end
    for y=3*L+3+1:4*L+4-1
        AM[y,y+1] = 1
    end
    for y=4*L+4+1:5*L+4-1
        AM[y,y+1] = 1
    end
    # the blue links
    for y=1:L
        AM[y,y+4] = 2
    end
    for y=L+1:2*L+1
        AM[y,y+5] = 2
    end
    for y=2*L+1+1:3*L+3-1
        AM[y,y+5] = 2
    end
    for y=3*L+3+1:4*L+4-1
        AM[y,y+4] = 2
    end
    # the green links
    for y=1:L
        AM[y,y+3] = 3
    end
    for y=L+1:2*L+1
        AM[y,y+4] = 3
    end
    for y=2*L+1+2:3*L+3
        AM[y,y+4] = 3
    end
    for y=3*L+3+2:4*L+4
        AM[y,y+3] = 3
    end
    return AM,numsites
end