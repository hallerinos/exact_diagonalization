using CSV, DelimitedFiles

data = hcat([[mod(idx-1,24)+1,idy,real(cp)] for (idx,cps) in enumerate(chiralProductExp[:,1:24]) for (idy,cp) in enumerate(cps)]...)
writedlm("./chiralproduct.csv",transpose(data),",")

data = hcat(tupls...)
writedlm("./chiralproduct_tuples.csv",transpose(data),",")