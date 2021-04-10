using SparseArrays, LinearAlgebra, KrylovKit
using LightGraphs, GraphPlot, UnicodePlots
using Revise
using JLD
using Plots
using LaTeXStrings
using Statistics
using ColorSchemes

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")

@time AM,numsites = buildAdjacencyMatrixSnowflakePBC()
@time SO, magZ = buildLatticeOperators(numsites, 0.5)

# SISJ = Array{SparseMatrixCSC}(undef,(3,numsites,numsites))
# for r=1:numsites, rpr=r+1:numsites, s=1:3
#     @info r,rpr,s
#     SISJ[s,r,rpr] = SO[s,r]*SO[s,rpr]
# end

# loop over all triangles (same orientation)
chiralOP = Vector(undef,48)
tupls = Vector(undef,48)
# r1 = 1; r2 = 2; r3 = 5
# chiralOP += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
# chiralOP += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
# chiralOP += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
# bulk down triangles
i = 1
for r1 in [1, 2, 13, 14, 15]
    r2 = r1+1
    r3 = r1+4
    chiralOP[i] = spzeros(SO[1,1].m,SO[1,1].n)
    chiralOP[i] += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
    chiralOP[i] += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
    chiralOP[i] += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
    tupls[i] = [r1, r2, r3]
    i+=1
end
for r1 in [4, 5, 6, 8, 9, 10, 11]
    r2 = r1+1
    r3 = r1+5
    chiralOP[i] = spzeros(SO[1,1].m,SO[1,1].n)
    chiralOP[i] += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
    chiralOP[i] += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
    chiralOP[i] += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
    tupls[i] = [r1, r2, r3]
    i+=1
end
# bulk up triangles
for r1 in [1, 2, 3, 14, 15]
    r2 = r1+4
    r3 = r1+3
    chiralOP[i] = spzeros(SO[1,1].m,SO[1,1].n)
    chiralOP[i] += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
    chiralOP[i] += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
    chiralOP[i] += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
    tupls[i] = [r1, r2, r3]
    i+=1
end
for r1 in [4, 5, 6, 7, 9, 10, 11]
    r2 = r1+5
    r3 = r1+4
    chiralOP[i] = spzeros(SO[1,1].m,SO[1,1].n)
    chiralOP[i] += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
    chiralOP[i] += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
    chiralOP[i] += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
    tupls[i] = [r1, r2, r3]
    i+=1
end
# the boundary terms
for tupl in [(1,17,2),(17,18,2),(2,18,3),(18,19,3),(3,8,7),(8,13,7),(7,13,12),(13,17,12),(12,1,16),(16,1,4),(16,4,19),(19,4,8),(18,19,3),(18,3,2),(17,18,2),(17,2,1),(13,17,12),(7,13,12),(8,13,7),(3,8,7),(19,4,8),(16,4,19),(16,1,4),(12,1,16)]
    r1 = tupl[1]
    r2 = tupl[2]
    r3 = tupl[3]
    chiralOP[i] = spzeros(SO[1,1].m,SO[1,1].n)
    chiralOP[i] += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
    chiralOP[i] += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
    chiralOP[i] += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])
    tupls[i] = [r1, r2, r3]
    i+=1
end
# total chirality
chiralOPSum = spzeros(SO[1,1].m,SO[1,1].n)
for cop in chiralOP
    chiralOPSum += cop
end

B = range(0,1,length=101)
ϵGS = Vector{Array}(undef, length(B))
magGS = Vector{Array}(undef, length(B))
sISJExp = Array{ComplexF64}(undef, (length(B),3,numsites,numsites))
chiralProductExp = Array{Array}(undef, (length(B), length(chiralOP)))
ϕGS = Vector{Array}(undef, length(B))
data = Nothing
# loop over all magnetic fields and compute observables
for (idx,bb) in enumerate(B)
    @info bb
    data = load("./saves/magField"*string(bb)*".jld")
    ham = data["ham"]

    λ = data["λ"]
    λ_ = [round(l; digits=5) for l in λ]
    λpr = unique(λ_)
    λₙ = [(l,count(x->x==l,λ_)) for l in λpr]

    # take the ground state manifold and ...
    ϕGS = hcat(data["ϕ"][1:λₙ[1][2]]...)

    submat = ϕGS'*chiralOPSum*ϕGS
    vals, vecs = eigen(submat)
    ϕGS *= vecs
    # ... switch to the eigenbasis of the chiral product
    for (idy,cop) in enumerate(chiralOP)
        # chiral product
        chiralProductExp[idx,idy] = [ϕGS[:,ii]'*cop*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
    end

    ϵGS[idx] = [ϕGS[:,ii]'*ham*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
    magGS[idx] = [ϕGS[:,ii]'*magZ*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]

    # # take the exp val of the correlators
    # for r=1:numsites, rpr=r+1:numsites, s=1:3
    #     @info bb,r,rpr
    #     sISJExp[idx,s,r,rpr] = ϕGS'*SISJ[s,r,rpr]*ϕGS
    # end
end