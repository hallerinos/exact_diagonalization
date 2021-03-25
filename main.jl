using LightGraphs
using GraphPlot
using UnicodePlots
using Arpack
using LuxurySparse
using Revise

includet("buildLatticeOperators.jl")
includet("buildAdjacencyMatrix.jl")
includet("buildHamiltonian.jl")
includet("diagonalize.jl")

spin = 0.5  # which local spin operator to take
L = [5, 4]  # defines the Lx, Ly cell
do_project = true  # project onto U(1) subspace, given by ...
projZ = -prod(L)/2+10  # the projection of the magnetization along the z-axis
# x, y, z couplings of the Heisenberg Hamiltonian
J₁ = -1.
J₂ = -1.
J₃ = -1.
# how many eigenvalues to take
nev = 20
# which tolerance in the eigenvalue problem
tol = 1e-12

@time SO = buildLatticeOperators(L, spin)
@time AM = buildAdjacencyMatrix(L)
@time ham = buildHamiltonian(SO, AM, J₁, J₂, J₃)
@time λ, ϕ = diagonalize(ham, nev, tol)
0;