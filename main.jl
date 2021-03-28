using SparseArrays, LinearAlgebra, KrylovKit
using LightGraphs, GraphPlot, UnicodePlots
using Revise
using JLD

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")
includet("buildHamiltonian.jl")
includet("diagonalize.jl")

spin = 0.5  # which local spin operator to take
L = [2, 3]  # defines the Lx, Ly cell -- will be ignored in case of snowflake configuration
do_project = false  # project onto U(1) subspace, given by ...
projZ = -prod(L)/2+10  # the projection of the magnetization along the z-axis -- is ignored for snowflake configuration
# the anisotropy term
Kᵤ = 0.0
# the DMI strength (use 1/sqrt(3) to make the DMI vector of unit norm)
𝐃 = 1.0
# x, y, z couplings of the Heisenberg Hamiltonian
J₁ = -0.5
J₂ = -0.5
J₃ = -0.5
# how many eigenvalues to take
nev = 20
# which tolerance in the eigenvalue problem
tol = 1e-12

@time AM,numsites = buildAdjacencyMatrixSnowflakePBC()  # build the adjacency matrix which defines the lattice connections
# @time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
display(gplot(DiGraph(AM), nodelabel=1:prod(numsites)))  # visualize the resulting graph
@time SO, magZ = buildLatticeOperators(numsites, spin, do_project=do_project, projZ=projZ)  # construct the lattice spin operators

# the magnetic field
# B = range(0,1,length=101)
B = range(0,1,length=101)
ii = 0
for bb in B
    @time ham = buildHamiltonian(SO, AM, J₁, J₂, J₃, Kᵤ, 𝐃, bb)    # construct the hamiltonian
    # # show(spy(real(ham)))
    # # println()
    @time λ, ϕ = diagonalize(ham, nev, tol)

    jldopen("system"*string(bb)*".jld", "w") do file
        write(file, "λ", λ)
        write(file, "ϕ", ϕ)
        write(file, "ham", ham)
    end

    ii+=1
end

0;  # just used to suppress the REPL output]