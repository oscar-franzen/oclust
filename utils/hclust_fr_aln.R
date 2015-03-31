args <- commandArgs(trailingOnly = TRUE)
args.working_dir <- args[1]

.libPaths(paste0(args.working_dir,"/bin/R.lib"))
library(seqinr)

print("hej")