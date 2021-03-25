include("spinOperators.jl")

function buildLatticeOperators(L::Array{Int64,1}, spin::Float64; do_project::Bool=false, projZ::Float64)
    d = Int64(2*spin+1)
    size = prod(L)
    D = d^size
    Sx, Sy, Sz, Id, Sm, Sp = spinOperators(spin)

    # create the sparse operators of the lattice
    SO = Array{SparseMatrixCSC}(undef,(3,size))
    magZ = spzeros(D,D)
    for r = 1:size
        # @info r
        IDL = sparse(I, d^(r-1), d^(r-1))
        IDR = sparse(I, d^(size-r), d^(size-r))
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
        for r = 1:size, i=1:3
            # @info r,i
            SO[i,r] = SO[i,r][idx,idx]
        end
    end

    return SO
end