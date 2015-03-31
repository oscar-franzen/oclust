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

z<-as.matrix(d)

# Workaround for NaN
z[is.na(z)]=1
d<-as.dist(z)

hc<-hclust(d,method=args.algo)

for (dist in c(0.005,1:10/100)) {
    cl<-cutree(hc,h=dist)

    df<-as.data.frame(cl)
    df<-data.frame(read_ID=rownames(df),cluster=df$cl)

    write.table(df, file=paste0(args.output_dir,"/MSA.",args.algo,".",dist,".hclust"), col.names=T, row.names=F, quote=F)
}