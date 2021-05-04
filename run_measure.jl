using SparseArrays, LinearAlgebra, KrylovKit
using LightGraphs, GraphPlot, UnicodePlots
using Revise
using JLD

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")
includet("buildHamiltonian.jl")
includet("diagonalize.jl")

spin = 0.5  # which local spin operator to take
L = [4, 4]  # defines the Lx, Ly cell -- will be ignored in case of snowflake configuration
do_project = false  # project onto U(1) subspace, given by ...
projZ = -prod(L)/2+10  # the projection of the magnetization along the z-axis -- is ignored for snowflake configuration
# the anisotropy term
Kᵤ = 0.0
# the DMI strength (use 1/sqrt(3) to make the DMI vector of unit norm)
𝐃 = 1.0
# x, y, z couplings of the Heisenberg Hamiltonian
J₁ = -0.5*abs(𝐃)
J₂ = -0.5*abs(𝐃)
J₃ = -0.5*abs(𝐃)
# J₁ = 1.0
# J₂ = 1.0
# J₃ = 1.0
# how many eigenvalues to take
nev = 20
# which tolerance in the eigenvalue problem
tol = 0.0

@time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
# @time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
# display(gplot(DiGraph(AM), nodelabel=1:prod(numsites)))  # visualize the resulting graph
@time SO, magZ = buildLatticeOperators(numsites, spin, do_project=do_project, projZ=projZ)  # construct the lattice spin operators

if !(@isdefined SISJ)
    SISJ = Array{SparseMatrixCSC}(undef,(3,3,numsites,numsites))
    SISIP1SJ = Array{SparseMatrixCSC}(undef,(3,3,3,numsites,numsites))
    SISJM1SJ = Array{SparseMatrixCSC}(undef,(3,3,3,numsites,numsites))
    for r=1:numsites, rpr=r+1:numsites, s=1:3, spr=1:3
        @info "create observables" r rpr s spr
        SISJ[s,spr,r,rpr] = SO[s,r]*SO[spr,rpr]
        for sprpr=1:3
            if r+1<=numsites
                SISIP1SJ[s,spr,sprpr,r,rpr] = SO[s,r]*SO[spr,r+1]*SO[sprpr,rpr]
            else
                SISIP1SJ[s,spr,sprpr,r,rpr] = spzeros(SO[1,1].m,SO[1,1].n)
            end
            if rpr==1
                SISJM1SJ[s,spr,sprpr,r,rpr] = spzeros(SO[1,1].m,SO[1,1].n)
            else
                SISJM1SJ[s,spr,sprpr,r,rpr] = SO[s,r]*SO[spr,rpr-1]*SO[sprpr,rpr]
            end
        end
    end
end

# the magnetic field
B = 20.0/39
# B = 0.
# B = range(0,0.5,length=2)
energies = []
λ = Nothing
ϕ = Nothing
ham = Nothing
@time ham = buildHamiltonian(SO, AM, J₁, J₂, J₃, Kᵤ, 𝐃, -B)    # construct the hamiltonian
# # show(spy(real(ham)))
# # println()

@time λ, ϕ = diagonalize(ham, nev, tol)
λ_ = [round(l; digits=5) for l in λ]
λpr = unique(λ_)
λₙ = [(l,count(x->x==l,λ_)) for l in λpr]

# ϕGS = hcat(ϕ[1:λₙ[1][2]]...)  # take the GS manifold
ϕGS = ϕ[1]  # to compare to MPS, take only the first state

ϵGS = Nothing
magGS = Nothing
SIExp = Array{Any}(undef, (3,numsites))
SISJExp = Array{Any}(undef, (3,3,numsites,numsites))
SISIP1SJExp = Array{Any}(undef, (3,3,3,numsites,numsites))
SISJM1SJExp = Array{Any}(undef, (3,3,3,numsites,numsites))

ϵGS = [ϕGS[:,ii]'*ham*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]
magGS =[ϕGS[:,ii]'*magZ*ϕGS[:,ii] for ii in 1:size(ϕGS,2)]

# take the exp val of the correlators
for r=1:numsites, s=1:3
    @info "measuring expectation values..." r s
    SIExp[s,r] = ϕGS'*SO[s,r]*ϕGS
    for rpr=r+1:numsites, spr=1:3
        # @info "measuring correlations..." r rpr s spr
        SISJExp[s,spr,r,rpr] = ϕGS'*SISJ[s,spr,r,rpr]*ϕGS
        for sprpr=1:3
            SISIP1SJExp[s,spr,sprpr,r,rpr] = ϕGS'*SISIP1SJ[s,spr,sprpr,r,rpr]*ϕGS
            SISJM1SJExp[s,spr,sprpr,r,rpr] = ϕGS'*SISJM1SJ[s,spr,sprpr,r,rpr]*ϕGS
        end
    end
end
SIExp[1,:]

open("local"*string(B)*".dat", "w") do io
    write(io, "site Re(X) Re(Y) Re(Z)\n")
    for r=1:numsites
        write(io, string(r), " ")
        for s=1:3
            write(io, string(real(SIExp[s,r])), " ");
        end
        write(io, "\n");
    end
end

open("corrB"*string(B)*".dat", "w") do io
    write(io, "site1 site2 Re(XX) Re(YY) Re(ZZ) Re(XY_Z) Re(X_YZ)\n")
    for r=1:numsites, rpr=r+1:numsites
        write(io, string(r), " ", string(rpr), " ")
        for s=1:3
            write(io, string(real(SISJExp[s,s,r,rpr])), " ");
        end
        write(io, string(real(SISIP1SJExp[1,2,3,r,rpr])), " ");
        write(io, string(real(SISJM1SJExp[1,2,3,r,rpr])), " ");
        write(io, "\n");
    end
end