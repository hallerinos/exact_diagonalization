using SparseArrays, LinearAlgebra, KrylovKit
using LightGraphs, GraphPlot, UnicodePlots
using Revise

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")
includet("buildHamiltonian.jl")
includet("diagonalize.jl")

spin = 0.5  # which local spin operator to take
L = [4,5]  # defines the Lx, Ly cell
do_project = true  # project onto U(1) subspace, given by ...
projZ = -prod(L)/2+10  # the projection of the magnetization along the z-axis
# x, y, z couplings of the Heisenberg Hamiltonian
J₁ = -1.
J₂ = -1.
J₃ = -1.
# the anisotropy term
Kᵤ = .1
# the DMI vector
D = 0.0
𝐃 = [D,D,D]
# the magnetic field
B = 0.
# how many eigenvalues to take
nev = 20
# which tolerance in the eigenvalue problem
tol = 1e-12

# @time AM,numsites = buildAdjacencyMatrixSnowflake(L)  # build the adjacency matrix which defines the lattice connections
@time AM,numsites = buildAdjacencyMatrixRegular(L)  # build the adjacency matrix which defines the lattice connections
# display(gplot(DiGraph(AM), nodelabel=1:prod(numsites)))  # visualize the resulting graph
@time SO = buildLatticeOperators(numsites, spin, do_project=do_project, projZ=projZ)  # construct the lattice spin operators
@time ham = buildHamiltonian(SO, AM, J₁, J₂, J₃, Kᵤ, 𝐃, B)    # construct the hamiltonian
# show(spy(real(ham)))
# println()
@time λ, ϕ = diagonalize(ham, nev, tol)

0;  # just used to suppress the REPL output