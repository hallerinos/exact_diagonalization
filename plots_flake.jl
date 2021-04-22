using Plots, Plots.PlotMeasures, LaTeXStrings
qxs = range(-π,π,length=101)
qys = range(-π,π,length=101)

for (idb,bb) in enumerate(B)
    plot_zs = Array{Float64}(undef, length(qxs), length(qys))
    for (idx,qx) in enumerate(qxs), (idy,qy) in enumerate(qys)
        @info idx,idy,qx,qy
        z = 0
        for r=1:numsites-1, rpr=r+1:numsites
            for s in sISJExp[idb,3,r,rpr]
                z += exp(1im*((xs[rpr]-xs[r])*qx+(ys[r]-ys[rpr])*qy))*s
            end
        end
        z = 2*real(z)/length(sISJExp[idb,3,1,2])^2
        plot_zs[idx,idy] = z
    end
    gr()
    p = heatmap(qxs,qys,plot_zs, c = :Purples_9, xlabel=L"q_x", ylabel=L"q_y", formatter = :plain, legend=false, tick_direction=:out, size=(200,200), dpi=300,xtickfontsize=8,ytickfontsize=8, left_margin=-12px, right_margin=-6px,top_margin=-6px,bottom_margin=-15px)
    # p = Plots.Image(plot_zs,(-π, π),(-π, π))
    fn = "structure_factor_PBC/B_"*string(bb)*".png"
    savefig(fn)
    # PGFPlots.Plots.save(fn, p)
    # display(p)
end