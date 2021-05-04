using SparseArrays, LinearAlgebra, KrylovKit, JLD

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")

L=[3,3]
@time AM,numsites = buildAdjacencyMatrixFlake(L,pbc=true)
@time SO, magZ = buildLatticeOperators(numsites, 0.5)
xs, ys = constructLatticeFlake(L)

SISJ = Array{SparseMatrixCSC}(undef,(3,numsites,numsites))
for r=1:numsites, rpr=1:numsites, s=1:3
    @info r,rpr,s
    SISJ[s,r,rpr] = SO[s,r]*SO[s,rpr]
end

# loop over all triangles (same orientation)
nStripe = L[1]
edge = L[1]
i = 1
chiralOP = Vector(undef,numsites)
# this loop only constructs the bulk chiral triple product
for k=1:2*L[2]-1
    for j=1:nStripe
        @info i
        op = spzeros(SO[1,1].m,SO[1,1].n)
        if k<L[2]
            if mod(i,edge)!=0
                # up triangles
                r1 = i
                r2 = i + 1
                r3 = i + nStripe + 1
                op += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
                op += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
                op += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
                # down triangles
                r1 = i
                r2 = i - nStripe + 1
                r3 = i + 1
                if r2 > 0
                    op += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
                    op += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
                    op += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
                end
            else
                nStripe += 1
                edge += nStripe
            end
        else
            if mod(i,edge)!=0
                # up triangles
                r1 = i
                r2 = i + 1
                r3 = i + nStripe + 1
                if r3<=numsites
                    op += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
                    op += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
                    op += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
                end
                # down triangles
                r1 = i
                r2 = i - nStripe + 1
                r3 = i + 1
                if r3 <= numsites
                    op += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
                    op += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
                    op += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
                end
            else
                nStripe -= 1
                edge += nStripe
            end
        end
        chiralOP[i] = op
        i+=1
    end
end
# total chirality of the bulk
chiralOPSum = spzeros(SO[1,1].m,SO[1,1].n)
for cop in chiralOP
    chiralOPSum += cop
end

B = range(0,1,length=11)
ϵGS = Vector{Array}(undef, length(B))
magGS = Vector{Array}(undef, length(B))
sISJExp = Array{Array}(undef, (length(B), size(SO,1), size(SO,2), size(SO,2)))
chiralProductExp = Array{Array}(undef, (length(B), length(chiralOP)))
chiralProductSumExp = Vector{Array}(undef, length(B))
ϕGS = Vector{Array}(undef, length(B))
data = Nothing
# loop over all magnetic fields and compute observables
for (idx,bb) in enumerate(B)
    @info bb
    data = load("./saves_flake_PBC/magField"*string(bb)*".jld")
    ham = data["ham"]

    λ = data["λ"]
    λ_ = [round(l; digits=5) for l in λ]
    λpr = unique(λ_)
    λₙ = [(l,count(x->x==l,λ_)) for l in λpr]

    ϕGS = hcat(data["ϕ"][1:λₙ[1][2]]...)  # take the ground state manifold and ...

    submat = ϕGS'*chiralOPSum*ϕGS
    vals, vecs = eigen(submat)
    ϕGS *= vecs  # ... switch to the eigenbasis of the chiral product operator
    for (idy,cop) in enumerate(chiralOP)
        chiralProductExp[idx,idy] = [ϕGS[:,ii]'*cop*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
    end

    chiralProductSumExp[idx] = [ϕGS[:,ii]'*chiralOPSum*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]

    ϵGS[idx] = [ϕGS[:,ii]'*ham*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
    magGS[idx] = [ϕGS[:,ii]'*magZ/numsites*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]

    # take the exp val of the correlators
    for r=1:numsites, rpr=1:numsites
        @info idx,bb,r,rpr
        for s=1:3
            submat = ϕGS'*SISJ[s,r,rpr]*ϕGS
            if size(submat,1)>1
                vals, vecs = eigen(submat)
            else 
                vals = submat
            end
            sISJExp[idx,s,r,rpr] = vals
        end
    end
end