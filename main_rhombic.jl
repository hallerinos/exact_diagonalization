using SparseArrays, LinearAlgebra, KrylovKit
using LightGraphs, GraphPlot, UnicodePlots
using Revise
using JLD

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")
includet("buildHamiltonian.jl")
includet("diagonalize.jl")

spin = 0.5  # which local spin operator to take
L = [4, 5]  # defines the Lx, Ly cell -- will be ignored in case of snowflake configuration
do_project = false  # project onto U(1) subspace, given by ...
projZ = -prod(L)/2+10  # the projection of the magnetization along the z-axis -- is ignored for snowflake configuration
# the anisotropy term
Kᵤ = 0.0
# the DMI strength (use 1/sqrt(3) to make the DMI vector of unit norm)
𝐃 = 1.0
# x, y, z couplings of the Heisenberg Hamiltonian
J₁ = -0.5*𝐃
J₂ = -0.5*𝐃
J₃ = -0.5*𝐃
# J₁ = 1.0
# J₂ = 1.0
# J₃ = 1.0
# how many eigenvalues to take
nev = 1
# which tolerance in the eigenvalue problem
tol = 0.0

@time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
# @time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
display(gplot(DiGraph(AM), nodelabel=1:prod(numsites)))  # visualize the resulting graph
@time SO, magZ = buildLatticeOperators(numsites, spin, do_project=do_project, projZ=projZ)  # construct the lattice spin operators

# the magnetic field
B = range(0,1,length=40)
# B = range(0,0.5,length=2)
energies = []
for bb in B
    @info bb
    λ = Nothing
    ϕ = Nothing
    ham = Nothing
    @time ham = buildHamiltonian(SO, AM, J₁, J₂, J₃, Kᵤ, 𝐃, -bb)    # construct the hamiltonian
    # # show(spy(real(ham)))
    # # println()
    @time λ, ϕ = diagonalize(ham, nev, tol)
    println(λ[1])

    jldopen("saves_rhombic_OBC/magField"*string(bb)*".jld", "w") do file
        write(file, "λ", λ)
        write(file, "ϕ", ϕ)
        write(file, "ham", ham)
    end
end

0;  # just used to suppress the REPL output 