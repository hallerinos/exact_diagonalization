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

SISJ = Array{SparseMatrixCSC}(undef,(3,numsites,numsites))
for r=1:numsites, rpr=1:numsites, s=1:3
    SISJ[s,r,rpr] = SO[s,r]*SO[s,rpr]
end

B = range(0,1,length=101)[1:2]
ϵGS = Vector{ComplexF64}(undef, length(B))
magGS = Vector{ComplexF64}(undef, length(B))
SISJexp = Array{ComplexF64}(undef, (length(B),3,numsites,numsites))
ϕGS = Vector{Array}(undef, length(B))
data = Nothing
# loop over all magnetic fields and compute observables
for (idx,bb) in enumerate(B)
    data = load("./saves/magField"*string(bb)*".jld")
    ham = data["ham"]
    ϕGS = data["ϕ"][1]

    ϵGS[idx] = ϕGS'*ham*ϕGS
    magGS[idx] = ϕGS'*magZ*ϕGS

    # take the exp val of the SISJ correlator
    for r=1:numsites, rpr=1:numsites, s=1:3
        SISJexp[idx,s,r,rpr] = ϕGS'*SISJ[s,r,rpr]*ϕGS
    end

    # take the chiral triple product exp val
    for r=1:numsites, rpr=1:numsites, s=1:3
        SISJexp[idx,s,r,rpr] = ϕGS'*SISJ[s,r,rpr]*ϕGS
    end
end

# energy plot
scatter(B,real(ϵGS),linecolor=:black,lw=1)
scatter!(ylabel=L"E_0")
scatter!(xlabel=L"B")
scatter!(legend=false)
scatter!(yrange=[-17,-11])
scatter!(xrange=[0,1])
fn = "energy.pdf"
savefig(fn)

# magnetization scatter
scatter(B,real(-magGS)/numsites,linecolor=:black,lw=1)
scatter!(ylabel=L"\langle M_z\rangle")
scatter!(xlabel=L"B")
scatter!(legend=false)
scatter!(yrange=[0,0.5])
scatter!(xrange=[0,1])
fn = "magnetization.pdf"
savefig(fn)