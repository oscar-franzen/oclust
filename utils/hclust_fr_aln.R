args <- commandArgs(trailingOnly = TRUE)
args.working_dir <- args[1]
args.input_file <- args[2]
args.output_dir <- args[3]
args.algo <- args[4]

.libPaths(paste0(args.working_dir,"/bin/R.lib"))

library(seqinr)

a<-read.alignment(args.input_file,format="fasta")
d<-dist.alignment(a,matrix="identity")

# Because the returned distance is the square root of the distance
m<-as.matrix(d)^2

d<-as.dist(m)
hc<-hclust(d,method="complete")

for (dist in 1:10/100) {
    cl<-cutree(hc,h=dist)
    write.table(as.data.frame(cl), file=paste0("./infernal.hclust.out/",file,".",dist,".hclust"), row.names=T, quote=F)
}