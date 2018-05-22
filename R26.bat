x <- read.csv("GSE59525_RPKM_Dqua.txt.gz",sep="\t")
load("x1_Dqua3CW")
S <- split(1:dim(x1)[1],x1[,1])
SS <- list(rep(NA,length(S)))
for (i in c(1:length(S)))
{
cat(i, " ")
SS[[i]] <- split(S[[i]],x1[S[[i]],3])
}
gtf <- read.csv("DQUA.v02.gff3",sep="\t",comment.char="#",header=F) #Dqua
gtf1 <- gtf[grep("gene_name",gtf[,9]),] #Dqua

methyl <- NULL
for (i in c(1:dim(x)[1]))
{
cat(i," ")
if (x[i,5]=="-") {is<-1} else {is<-2} #Dqua
jj <- grep(x[i,1],gtf1[,9])[1] #Dqua
ii <- match(gtf1[jj,1],names(S)) #Dqua
index<-  SS[[ii]][[is]][x1[SS[[ii]][[is]],2]>=gtf1[jj,4] & x1[SS[[ii]][[is]],2]<=gtf1[jj,5]] #Dqua
methyl <- rbind(methyl,data.frame(x[i,1],t(colSums(x1[index,4:5]))))
}
save(file="methyl_Dqua3CW",methyl)
