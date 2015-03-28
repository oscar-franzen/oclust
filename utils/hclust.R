args <- commandArgs(trailingOnly = TRUE)

arg.target_matrix<-args[1]
arg.output_file<-args[2]
arg.method<-args[3]

f<-read.table(arg.target_matrix,header=T)

m<-as.matrix(f)
d<-as.dist(m)

hc<-hclust(d,method=arg.method)
hc.c<-cutree(hc,h=0.01)

df<-data.frame(read=rownames(f),cluster=hc.c)

write.table(df, file=paste0(arg.output_file,".0.01.hclust"), row.names=F, quote=F)

hc<-hclust(d,method=arg.method)
hc.c<-cutree(hc,h=0.02)

df<-data.frame(read=rownames(f),cluster=hc.c)

write.table(df, file=paste0(arg.output_file,".0.02.hclust"), row.names=F, quote=F)

hc<-hclust(d,method=arg.method)
hc.c<-cutree(hc,h=0.03)

df<-data.frame(read=rownames(f),cluster=hc.c)

write.table(df, file=paste0(arg.output_file,".0.03.hclust"), row.names=F, quote=F)

hc<-hclust(d,method=arg.method)
hc.c<-cutree(hc,h=0.04)

df<-data.frame(read=rownames(f),cluster=hc.c)

write.table(df, file=paste0(arg.output_file,".0.04.hclust"), row.names=F, quote=F)
