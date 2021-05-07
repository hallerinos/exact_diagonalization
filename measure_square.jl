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

@time AM, numsites = buildAdjacencyMatrixRegular([4,4])
@time SO, magZ = buildLatticeOperators(numsites, 0.5)

SISJ = Array{SparseMatrixCSC}(undef,(3,3,numsites,numsites))
for r=1:numsites, rpr=r+1:numsites, s=1:3, spr=1:3
    @info "create SISJ" r rpr s spr
    SISJ[s,spr,r,rpr] = SO[s,r]*SO[spr,rpr]
end

SISJSK = Array{SparseMatrixCSC}(undef,(3,3,3,numsites,numsites,numsites))
for r=1:numsites, rpr=r+1:numsites, rprpr=rpr+1:numsites, s=1:3, spr=1:3, sprpr=1:3
    @info "create SISJSK" r rpr s spr sprpr
    SISJSK[s,spr,sprpr,r,rpr,rprpr] = SO[s,r]*SO[spr,rpr]*SO[sprpr,rprpr]
end
## 
B = range(0,1,length=40)
ϵGS = Vector{Array}(undef, length(B))
magGS = Vector{Array}(undef, length(B))
sIExp = Array{Array}(undef, (length(B),3,numsites))
sISJExp = Array{Array}(undef, (length(B),3,3,numsites,numsites))
sISJSKExp = Array{Array}(undef, (length(B),3,3,3,numsites,numsites,numsites))
ϕGS = Nothing
data = Nothing
# loop over all magnetic fields and compute observables
for (idx,bb) in enumerate(B)
    @info "Measure expectation values for " bb
    data = load("./saves_square_OBC/magField"*string(bb)*".jld")
    ham = data["ham"]

    λ = data["λ"]
    λ_ = [round(l; digits=5) for l in λ]
    λpr = unique(λ_)
    λₙ = [(l,count(x->x==l,λ_)) for l in λpr]

    # take the ground state manifold and ...
    ϕGS = hcat(data["ϕ"][1:λₙ[1][2]]...)

    ϵGS[idx] = [ϕGS[:,ii]'*ham*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
    magGS[idx] = [ϕGS[:,ii]'*magZ*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]

    # take the exp val of the correlators
    for r=1:numsites, s=1:3
        # @info "measuring local expectation values..." r s
        sIExp[idx,s,r] = ϕGS'*SO[s,r]*ϕGS
    end

    path = "saves_square_OBC/observables/B"*string(bb)
    if !isdir(path)
        mkpath(path)
    end
    open(path*"/local.dat", "w") do io
        write(io, "i <Sx_i> <Sy_i> <Sz_i>\n")
        for i=1:numsites
            write(io, string(i))
            for s=1:3
                write(io, " ", string(real(sIExp[idx,s,i][1])));
            end
            write(io, "\n");
        end
    end

    for r=1:numsites, rpr=r+1:numsites, s=1:3, spr=1:3
        sISJExp[idx,s,spr,r,rpr] = ϕGS'*SISJ[s,spr,r,rpr]*ϕGS
        for rprpr=rpr+1:numsites,sprpr=1:3
            # @info "measuring correlations..." r rpr s spr sprpr
            sISJSKExp[idx,s,spr,sprpr,r,rpr,rprpr] = ϕGS'*SISJSK[s,spr,sprpr,r,rpr,rprpr]*ϕGS
        end
    end
    
    open(path*"/corr2.dat", "w") do io
        write(io, "i j <Sx_iSx_j> <Sx_iSy_j> <Sx_iSz_j> <Sy_iSy_j> <Sy_iSz_j> <Sz_iSz_j>\n")
        for i=1:numsites, j=i+1:numsites
            write(io, string(i), " ", string(j))
            write(io, " ", string(real(sISJExp[idx,1,1,i,j][1])));
            write(io, " ", string(real(sISJExp[idx,1,2,i,j][1])));
            write(io, " ", string(real(sISJExp[idx,1,3,i,j][1])));
            write(io, " ", string(real(sISJExp[idx,2,2,i,j][1])));
            write(io, " ", string(real(sISJExp[idx,2,3,i,j][1])));
            write(io, " ", string(real(sISJExp[idx,3,3,i,j][1])));
            write(io, "\n");
        end
    end

    open(path*"/corr3.dat", "w") do io
        write(io, "i j k <Sx_iSy_jSz_k> <Sz_iSx_jSy_k> <Sy_iSz_jSx_k>\n")
        for i=1:numsites-1, j=i+1:numsites-1, k=j+1:numsites
            write(io, string(i), " ", string(j), " ", string(k))
            write(io, " ", string(real(sISJSKExp[idx,1,2,3,i,j,k][1])));
            write(io, " ", string(real(sISJSKExp[idx,3,1,2,i,j,k][1])));
            write(io, " ", string(real(sISJSKExp[idx,2,3,1,i,j,k][1])));
            write(io, "\n");
        end
    end
end