using LightGraphs
using GraphPlot
using UnicodePlots
using Arpack
using LuxurySparse
using Revise

includet("buildLatticeOperators.jl")

spin = 0.5
L = [4, 4]
do_project = false
projZ = prod(L)/2-5
J₁ = -1.
J₂ = -1.
J₃ = -1.
nev = 20
tol = 1e-12

@time buildLatticeOperators(L, spin, nev, tol, J₁, J₂, J₃, do_project=do_project, projZ=projZ)

0;