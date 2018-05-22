x <- read.csv("GSE59525_RPKM_Pcan.txt.gz",sep="\t")
load("x1_Pcan_S")
S <- split(1:dim(x1)[1],x1[,1])
SS <- list(rep(NA,length(S)))
for (i in c(1:length(S)))
{
cat(i, " ")
SS[[i]] <- split(S[[i]],x1[S[[i]],3])
}
methyl <- NULL
for (i in c(1:dim(x)[1]))
{
cat(i," ")
if (x[i,5]=="+") {is<-1} else {is<-2}
index<-  SS[[x[i,2]]][[is]][x1[SS[[x[i,2]]][[is]],2]>=x[i,3] & x1[SS[[x[i,2]]][[is]],2]<=x[i,4]]
methyl <- rbind(methyl,data.frame(x[i,1],t(colSums(x1[index,4:5]))))
}
save(file="methyl_Pcan_S",methyl)
