function buildHamiltonian(SO::Array{SparseMatrixCSC,2}, AM::Array{Float64,2}, J₁::Float64, J₂::Float64, J₃::Float64, Kᵤ::Float64, 𝐃::Array{Float64,1}, B::Float64)
    # build the hamiltonian matrix
    ham = spzeros(SO[1,1].m,SO[1,1].n)
    for r=1:size(SO,2), rpr=1:size(SO,2)
        # println(r, " ", rpr)
        # if adjacency matrix has a nonzero entry, allow the interactions
        if AM[r, rpr] != 0
            # the Heisenberg terms
            ham += J₁*SO[1,r]*SO[1,rpr]
            ham += J₂*SO[2,r]*SO[2,rpr]
            ham += (J₃ + Kᵤ)*SO[3,r]*SO[3,rpr]
            # the magnetic field
            ham += B*SO[3,r]
        end
        if AM[r, rpr] == 1
            # the red DM terms
            ham += 𝐃[1] * (SO[1:3,r] × SO[1:3,rpr])[2]
        end
        if AM[r, rpr] == 2
            # the red DM terms
            ham += 0.5 * 𝐃[2] * (SO[1:3,r] × SO[1:3,rpr])[2]
            ham += sqrt(3)/2 * 𝐃[2] * (SO[1:3,r] × SO[1:3,rpr])[1]
        end
        if AM[r, rpr] == 3
            # the red DM terms
            ham += -0.5 * 𝐃[2] * (SO[1:3,r] × SO[1:3,rpr])[2]
            ham += sqrt(3)/2 * 𝐃[2] * (SO[1:3,r] × SO[1:3,rpr])[1]
        end
    end
    return ham
end