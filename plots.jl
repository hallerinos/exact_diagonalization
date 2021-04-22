# data = load("./saves/magField"*string(0.12)*".jld")
# ham = data["ham"]

# λ = data["λ"]
# λ_ = [round(l; digits=5) for l in λ]
# λpr = unique(λ_)
# λₙ = [(l,count(x->x==l,λ_)) for l in λpr]

# # take the ground state manifold and ...
# ϕGS = hcat(data["ϕ"][1:λₙ[1][2]]...)
# # ... switch to the eigenbasis of the chiral product
# # submat = ϕGS'*chiralOP*ϕGS
# # vals, vecs = eigen(submat)
# # ϕGS *= vecs

# chiralOP = spzeros(SO[1,1].m,SO[1,1].n)
# r1 = 1; r2 = 2; r3 = 5;
# r1 = 2; r2 = 3; r3 = 6;
# chiralOP += SO[1,r1]*(SO[2,r2]*SO[3,r3] - SO[3,r2]*SO[2,r3])
# chiralOP += SO[2,r1]*(SO[3,r2]*SO[1,r3] - SO[1,r2]*SO[3,r3])
# chiralOP += SO[3,r1]*(SO[1,r2]*SO[2,r3] - SO[2,r2]*SO[1,r3])

# ϕGS[:,1]'*chiralOP*ϕGS[:,1]*38/pi


# # energy plot
# scatter(B,real(ϵGS),linecolor=:black,lw=1)
# scatter!(ylabel=L"E_0")
# scatter!(xlabel=L"B")
# scatter!(legend=false)
# scatter!(yrange=[-17,-11])
# scatter!(xrange=[0,1])
# fn = "energy.pdf"
# savefig(fn)

# # magnetization scatter plot
# scatter(B,-real(magGS)/numsites,lw=1, label=L"\langle M_z\rangle")
# scatter!(B,-real(chiralProductExp)*38/pi,lw=1, label=L"\langle Q_\psi\rangle")
# scatter!(xlabel=L"B")
# scatter!(yrange=[0,0.5])
# scatter!(xrange=[0,1])
# fn = "magnetization.pdf"
# savefig(fn)

colors = ["black","red","blue","green","yellow","gray"]

magData = hcat([[B[idx],-real(mag),colors[idy]] for (idx,mags) in enumerate(magGS) for (idy,mag) in enumerate(mags)]...)
scatter(magData[1,:],magData[2,:], color=magData[3,:], legend=false, xlabel=L"B", ylabel=L"M_z = \frac{1}{N}\sum_{i\in\{1,\dots,N\}}\langle{S_{z,i}}\rangle")
fn = "magnetization.pdf"
savefig(fn)

ϵData = hcat([[B[idx],real(ϵ),colors[idy]] for (idx,ϵs) in enumerate(ϵGS) for (idy,ϵ) in enumerate(ϵs)]...)
scatter(ϵData[1,:],ϵData[2,:], color=ϵData[3,:], legend=false, xlabel=L"B", ylabel=L"E_0")
fn = "energy.pdf"
savefig(fn)

cpData = hcat([[B[idx],real(cp),colors[idy]] for (idx,cps) in enumerate(chiralProductSumExp) for (idy,cp) in enumerate(cps)]...)
# scatter(cpData[1,:],-48*cpData[2,:]/pi, color=cpData[3,:], legend=false, xlabel=L"B", ylabel=L"\langle Q_\psi\rangle = 1/\pi\sum_{(i,j,k)\in\{\Delta,\nabla\}}\langle{S}_i\cdot({S}_j\times{S}_k)\rangle")
# fn = "Q_psi_diag.pdf"
# savefig(fn)
scatter(cpData[1,:],cpData[2,:]/pi, color=cpData[3,:], legend=false, xlabel=L"B", ylabel=L"\langle Q_\psi\rangle = 48/\pi\langle{S}_1\cdot({S}_2\times{S}_5)\rangle")
fn = "Q_psi_2_diag.pdf"
savefig(fn)

# scatter([(B[idx],-real(cp)/pi) for (idx,cps) in enumerate(chiralProductExp) for cp in cps], color=:black, legend=false, xlabel=L"B", ylabel=L"\langle Q_\psi\rangle = \frac{1}{\pi}\sum_{(i,j,k)\in\{\Delta,\nabla\}}\langle{S}_i\cdot({S}_j\times{S}_k)\rangle")
# fn = "Q_psi.pdf"
# savefig(fn)