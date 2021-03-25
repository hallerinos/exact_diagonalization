function buildHamiltonian(SO::Array{SparseMatrixCSC,2}, AM::Array{Float64,2}, J₁::Float64, J₂::Float64, J₃::Float64)
    # build the hamiltonian matrix
    # ham = spzeros(D,D)
    ham = spzeros(SO[1,1].m,SO[1,1].n)
    for r=1:size(SO,2), rpr=1:size(SO,2)
        # println(r, " ", rpr)
        # if adjacency matrix has a nonzero entry, allow the interactions
        if AM[r, rpr] != 0
            ham += J₁*SO[1,r]*SO[1,rpr]
            ham += J₂*SO[2,r]*SO[2,rpr]
            ham += J₃*SO[3,r]*SO[3,rpr]
        end
    end
    # show(spy(real(ham)))
    # println()
    return ham
end