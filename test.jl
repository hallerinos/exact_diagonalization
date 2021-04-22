using PGFPlots

p = Plots.Histogram(rand(10))
save("myfile.pdf", p)