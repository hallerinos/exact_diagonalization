function buildAdjacencyMatrixFlake(L::Array{Int64,1}; pbc::Bool=true)
    numsites = (L[2]-1)^2 + L[1]*(2*L[2]-1)

    # build the adjacency matrix
    AM = zeros(numsites, numsites)

    # build the red links
    # pbc terms
    if pbc==true
        # red terms
        i = 0
        idxs = []
        for j=0:L[2]-1
            i += L[1]+j
            append!(idxs,i)
        end
        for i=1:length(idxs)
            # @info idxs[i],size(AM,2)-idxs[end-i+1]
            AM[idxs[i],end-idxs[end-i+1]+1] = 1
        end
        idxs = []
        i = 1
        for j=0:L[2]-2
            append!(idxs,i)
            i += L[1]+j
        end
        for i=1:length(idxs)
            # @info size(AM,1)-idxs[i]+1,idxs[end-i+1]
            AM[end-idxs[i]+1,idxs[end-i+1]] = 1
        end

        # green terms
        idxs = []
        i=1
        for j=0:L[2]-1
            # @info i
            append!(idxs,i)
            i += L[1] + j
        end
        for i=1:length(idxs)
            # @info idxs[i], size(AM,2)-idxs[end-i+1]+1
            AM[end-idxs[i]+1,idxs[end-i+1]] = 2
        end
        idxs = []
        for j=2:L[1]
            append!(idxs,j)
        end
        for i=1:length(idxs)
            # @info idxs[i],size(AM,2)-idxs[end-i+1]
            AM[end-idxs[end-i+1]+1,idxs[i]] = 2
        end

        # blue terms
        idxs = []
        for j=1:L[1]
            append!(idxs,j)
        end
        for i=1:length(idxs)
            # @info idxs[i],size(AM,2)-idxs[end-i+1]
            AM[end-idxs[end-i+1]+1,idxs[i]] = 3
        end
        idxs = []
        i = L[1]
        for j=1:L[2]-1
            i += L[1]+j
            append!(idxs,i)
        end
        for i=1:length(idxs)
            AM[end-idxs[end-i+1]+1,idxs[i]] = 3
        end
    end

    # obc terms
    i = 1
    nStripe = L[1]
    edge = L[1]
    for k=1:L[2]-1
        for j=1:nStripe
            if mod(i,edge)!=0
                # red terms
                AM[i,i+1] = 1
                AM[end-i,end-i+1] = 1
                # green terms
                AM[i,i+nStripe+1] = 2
                AM[end-i-nStripe,end-i+1] = 2
                # blue terms
                AM[i,i+nStripe] = 3
                AM[end-i-nStripe+1,end-i+1] = 3
            else
                # green terms
                AM[i,i+nStripe+1] = 2
                AM[end-i-nStripe,end-i+1] = 2
                # blue terms
                AM[i,i+nStripe] = 3
                AM[end-i-nStripe+1,end-i+1] = 3
                nStripe += 1
                edge += nStripe
            end
            i += 1
        end
    end
    for j=i:i+L[1]+L[2]-3
        AM[j,j+1] = 1
    end

    return AM,numsites
end

function buildAdjacencyMatrixRegular(L::Array{Int64,1})
    numsites = prod(L)
    # build the adjacency matrix
    AM = zeros(numsites, numsites)
    for r=1:numsites
        # the red links
        mod(r,L[1])!=0 && r+1<=numsites ? AM[r, r+1] = 1 : Nothing
        # the green links
        mod(r,L[1])!=0 && r+L[1]+1<=numsites ? AM[r, r+L[1]+1] = 2 : Nothing
        # the blue links
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

function constructLatticeFlake(L::Array{Int64,1})
    δ₁ = [1.,0.]
    δ₂ = [0.5,-sqrt(3.)/2.]
    δ₃ = [-0.5,-sqrt(3.)/2.]

    xs = []
    ys = []
    pos = [0.,0.]
    nStripe = L[1]
    edge = L[1]
    tp = pos
    i = 1
    for k=1:2*L[2]-1
        for j=1:nStripe
            append!(xs,pos[1])
            append!(ys,pos[2])
            if k<L[2]
                if mod(i,edge)!=0
                    pos += δ₁
                else
                    pos = tp + δ₃
                    tp = pos
                    nStripe += 1
                    edge += nStripe
                end
            else
                if mod(i,edge)!=0
                    pos += δ₁
                else
                    pos = tp + δ₂
                    tp = pos
                    nStripe -= 1
                    edge += nStripe
                end
            end
            i+=1
        end
    end
    # p = plot()
    # for i=1:length(xs)
    #     for j=1:length(xs)
    #         if AM[i,j]==1
    #             quiver!(xs[i:i],ys[i:i],quiver=(xs[j:j]-xs[i:i],ys[j:j]-ys[i:i]),color=:red,linestyle=:dot)
    #         end
    #         if AM[i,j]==2
    #             quiver!(xs[i:i],ys[i:i],quiver=(xs[j:j]-xs[i:i],ys[j:j]-ys[i:i]),color=:green,linestyle=:dash)
    #         end
    #         if AM[i,j]==3
    #             quiver!(xs[i:i],ys[i:i],quiver=(xs[j:j]-xs[i:i],ys[j:j]-ys[i:i]),color=:blue,linestyle=:dashdot)
    #         end
    #     end
    # end
    # scatter!(xs,ys,legend=false)
    # display(p)

    return xs, ys
end