using Plots, Plots.PlotMeasures, LaTeXStrings
qxs = range(-π,π,length=101)
qys = range(-π,π,length=101)

for (idb,bb) in enumerate(B)
    for (ids,s) in enumerate(sISJExp[idb,2,1,2])
        zs = Array{Float64}(undef, length(qxs), length(qys))
        for (idx,qx) in enumerate(qxs), (idy,qy) in enumerate(qys)
            @info idx,idy,qx,qy
            z = 0
            for r=1:numsites-1, rpr=r+1:numsites
                z += exp(1im*((xs[rpr]-xs[r])*qx+(ys[r]-ys[rpr])*qy))*sISJExp[idb,2,r,rpr][ids]
            end
            z = 2*real(z)/numsites
            zs[idx,idy] = z
        end
        p = heatmap(qxs,qys,zs,c = :Purples_9, xlabel=L"q_x", ylabel=L"q_y", formatter = :plain, legend=true, tick_direction=:out, size=(200,200), dpi=300,xtickfontsize=8,ytickfontsize=8, left_margin=-12px, right_margin=8px,top_margin=0px,bottom_margin=-15px)
        fn = "structure_factor_YY_PBC/B_"*string(bb)*"ev_"*string(ids)*".png"
        savefig(fn)
    end
end

# cmax = maximum(plot_zs)
# cmin = minimum(plot_zs)
# absmax = maximum([cmax,cmin])
# cmax = cmax/absmax
# cmin = cmin/absmax
# for (idb,bb) in enumerate(B)
#     p = heatmap(qxs,qys,plot_zs[idb,:,:], c = :Purples_9, xlabel=L"q_x", ylabel=L"q_y", formatter = :plain, legend=true, tick_direction=:out, size=(200,200), dpi=300,xtickfontsize=8,ytickfontsize=8, left_margin=-12px, right_margin=8px,top_margin=0px,bottom_margin=-15px)
#     # p = heatmap(qxs,qys,abs.(plot_zs[idb,:,:])/absmax, c = :Purples_9, xlabel=L"q_x", ylabel=L"q_y", formatter = :plain, legend=true, tick_direction=:out, size=(200,200), dpi=300,xtickfontsize=8,ytickfontsize=8, left_margin=-12px, right_margin=8px,top_margin=0px,bottom_margin=-15px,clims=(0,cmax))
#     # p = Plots.Image(plot_zs,(-π, π),(-π, π))
#     fn = "structure_factor_PBC/B_"*string(bb)*".png"
#     savefig(fn)
# end