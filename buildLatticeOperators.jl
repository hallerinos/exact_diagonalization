using SparseArrays, LinearAlgebra, KrylovKit
include("spinOperators.jl")

function buildLatticeOperators(L::Array{Int64,1}, spin::Float64, nev::Int64, tol::Float64, J₁::Float64, J₂::Float64, J₃::Float64; do_project::Bool=false, projZ::Float64)
    d = Int64(2*spin+1)
    size = prod(L)
    D = d^size
    Sx, Sy, Sz, Id, Sm, Sp = spinOperators(spin)

    # create the sparse operators of the lattice
    SX = Vector{SparseMatrixCSC}(undef,size)
    SY = Vector{SparseMatrixCSC}(undef,size)
    SZ = Vector{SparseMatrixCSC}(undef,size)
    magZ = spzeros(D,D)
    for r = 1:size
        println(r)
        IDL = sparse(I, d^(r-1), d^(r-1))
        IDR = sparse(I, d^(size-r), d^(size-r))
        SX[r] = kron(IDL, Sx, IDR)
        SY[r] = kron(IDL, Sy, IDR)
        SZ[r] = kron(IDL, Sz, IDR)
        magZ += SZ[r]
    end
    ID = sparse(I, d^size, d^size)

    Dproj = D
    if do_project
        # write the total magnetization
        magZVec = diag(magZ)
        idx = findall(x->x==projZ, magZVec)
        Dproj = length(idx)
        # project all operators on the chosen subspace of the total magnetization
        ID = ID[idx,idx]
        for r = 1:size
            println(r)
            SX[r] = SX[r][idx,idx]
            SY[r] = SY[r][idx,idx]
            SZ[r] = SZ[r][idx,idx]
        end
    end

    # build the adjacency matrix
    AM = spzeros(size, size)
    for r=1:size
        mod(r,L[1])!=0 && r+1<=size ? AM[r, r+1] = 1 : Nothing
        r+L[1]<=size ? AM[r, r+L[1]] = 1 : Nothing
        mod(r,L[1])!=0 && r+L[1]+1<=size ? AM[r, r+L[1]+1] = 1 : Nothing
    end
    # show(spy(AM))
    # display(gplot(DiGraph(AM), nodelabel=1:size))

    # build the hamiltonian matrix
    # ham = spzeros(D,D)
    ham = spzeros(Dproj,Dproj)
    for r=1:size, rpr=1:size
        println(r, " ", rpr)
        # if adjacency matrix has a nonzero entry, allow the interactions
        if AM[r, rpr] != 0
            ham += J₁*SX[r]*SX[rpr]
            ham += J₂*SY[r]*SY[rpr]
            ham += J₃*SZ[r]*SZ[rpr]
        end
    end
    # show(spy(real(ham)))
    # println()

    nev_ = min(ham.m, nev)
    @time λ, ϕ = eigsolve(ham, nev, :SR, eltype(ham), issymmetric=true, krylovdim=max(nev_+10,4), tol=tol)
    # @time λ, ϕ = eigs(ham, nev=nev, which=:SR)
    # @time λ, ϕ = eigen(Matrix(ham))
    λ_ = [round(l; digits=5) for l in λ]
    λpr = unique(λ_)
    λₙ = [(l,count(x->x==l,λ_)) for l in λpr]
    print(λₙ)
    return λ, ϕ
end