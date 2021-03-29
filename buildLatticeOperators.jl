include("spinOperators.jl")

function buildLatticeOperators(numsites::Int64, spin::Float64; do_project::Bool=false, projZ::Float64=0.0)
    d = Int64(2*spin+1)
    D = d^numsites
    Sx, Sy, Sz, Id, Sm, Sp = spinOperators(spin)

    # create the sparse operators of the lattice
    SO = Array{SparseMatrixCSC}(undef,(3,numsites))
    magZ = spzeros(D,D)
    for r = 1:numsites
        # @info r
        IDL = sparse(I, d^(r-1), d^(r-1))
        IDR = sparse(I, d^(numsites-r), d^(numsites-r))
        SO[1,r] = kron(IDL, Sx, IDR)
        SO[2,r] = kron(IDL, Sy, IDR)
        SO[3,r] = kron(IDL, Sz, IDR)
        magZ += SO[3,r]
    end

    Dproj = D
    if do_project
        # write the total magnetization
        magZVec = diag(magZ)
        idx = findall(x->x==projZ, magZVec)
        Dproj = length(idx)
        # project all operators on the chosen subspace of the total magnetization
        for r = 1:numsites, i=1:3
            # @info r,i
            SO[i,r] = SO[i,r][idx,idx]
        end
    end

    return SO, magZ
end