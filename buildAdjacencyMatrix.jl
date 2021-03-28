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

function buildAdjacencyMatrixSnowflake()
    L = 3
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

function buildAdjacencyMatrixSnowflakePBC()
    L = 3
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
    AM[3,8] = 1
    AM[7,13] = 1
    AM[12,17] = 1
    AM[16,1] = 1
    AM[19,4] = 1
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
    AM[19,8] = 2
    AM[16,4] = 2
    AM[12,1] = 2
    AM[17,2] = 2
    AM[18,3] = 2
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
    AM[8,7] = 3
    AM[13,12] = 3
    AM[17,1] = 3
    AM[18,2] = 3
    AM[19,3] = 3
    return AM,numsites
end