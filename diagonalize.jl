function diagonalize(ham::SparseMatrixCSC{Complex{Float64},Int64}, nev::Int64, tol::Float64)
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