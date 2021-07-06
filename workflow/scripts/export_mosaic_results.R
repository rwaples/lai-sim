
library('MOSAIC')
model_results = snakemake@input[["model_results"]]
la_results = snakemake@input[["la_results"]]
mosaic_input_dir = snakemake@params[["input_dir"]]
output_path = snakemake@output[["path"]]

load(model_results)
load(la_results)
# get local ancestry at each SNP
local_pos=grid_to_pos(localanc, mosaic_input_dir, g.loc, chrnos)
# convert to array and then export
arr = array(unlist(local_pos, use.names=FALSE), c(3, 4, 20000))
save(arr, file = output_path)




# for each position and haplotype, find the ancestry with the largest probability
# threshold = 0.9
# apply threshold
# arr[arr<threshold] <- NA
# calls = apply(arr, c(2,3), function(x) which(x == max(x, na.rm=TRUE))[1])
# calls = t(calls)
# zerolen = apply(calls, c(1,2), is.na)
# calls[zerolen] <- 0
