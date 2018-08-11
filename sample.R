#Dqua methyl by pca
x1 <- read.csv("GSM1438857_Dqua_Synthetic_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua_S",x1)
x1 <- read.csv("GSM1438858_Dqua1AQ_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua1AQ",x1)
x1 <- read.csv("GSM1438859_Dqua2AQ_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua2AQ",x1)
x1 <- read.csv("GSM1438860_Dqua3AQ_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua3AQ",x1)
x1 <- read.csv("GSM1438861_Dqua1CW_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua1CW",x1)
x1 <- read.csv("GSM1438862_Dqua3CW_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua3CW",x1)
x1 <- read.csv("GSM1438863_Dqua3DW_CX.txt.gz",sep="\t",header=F)
save(file="x1_Dqua3DW",x1)

#Execute the following outside R
#nohup R --slave --vanilla < R21.bat > nohup21.out &
#nohup R --slave --vanilla < R22.bat > nohup22.out &
#nohup R --slave --vanilla < R23.bat > nohup23.out &
#nohup R --slave --vanilla < R24.bat > nohup24.out &
#nohup R --slave --vanilla < R25.bat > nohup25.out &
#nohup R --slave --vanilla < R26.bat > nohup26.out &
#nohup R --slave --vanilla < R27.bat > nohup27.out &


files<- list.files(pattern="methyl_D") 
files <-  files[c(7,2,5,6,1,3,4)]
methyl_all <-NULL
for (i in c(1:length(files)))
{
load(files[i])
methyl_all <- cbind(methyl_all,methyl[,2] / rowSums(methyl[,2:3]))
}
methyl_all[is.na(methyl_all)] <- 0

pca <- prcomp(scale(methyl_all))
P <- pchisq(scale(pca$x[,1])^2,1,lower.tail=F)
index <- p.adjust(P,"BH")<0.01
methyl[index,1]

#Dqua mRNA expression by pca
x <- read.csv("GSE59525_RPKM_Dqua.txt.gz",sep="\t")
pca1 <- prcomp(scale(x[,7:19]))
P1 <- pchisq(scale(pca1$x[,4])^2,1,lower.tail=F)
index1 <- p.adjust(P1,"BH")<0.01
x[index1,1]

#Dqua by HOSVD
require(rTensor)
Z <- array(NA,c(7,13,dim(x)[1]))
for (i in c(1:dim(x)[1]))
{
cat(i," ")
Z[,,i] <- outer(methyl_all[i,],data.matrix(x[i,7:19]),"*")
}
Z<-apply(Z,c(1,2),scale)
Z <- aperm(Z,c(2,3,1))
HOSVD <- hosvd(as.tensor(Z))
P3 <- pchisq(scale(HOSVD$U[[3]][,11])^2,1,lower.tail=F)
index3 <- p.adjust(P3,"BH")<0.01
methyl[index3,1]

#Pcan methyl by pca
x1 <- read.csv("GSM1438850_Pcan_Synthetic_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan_S",x1)
x1 <- read.csv("GSM1438851_Pcan21Q_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan21Q",x1)
x1 <- read.csv("GSM1438852_Pcan43Q_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan43Q",x1)
x1 <- read.csv("GSM1438853_Pcan75Q_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan75Q",x1)
x1 <- read.csv("GSM1438854_Pcan26W_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan26W",x1)
x1 <- read.csv("GSM1438855_Pcan42W_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan42W",x1)
x1 <- read.csv("GSM1438856_Pcan76W_CX.txt.gz",sep="\t",header=F)
save(file="x1_Pcan76W",x1)

#Execute the following outside R
#nohup R --slave --vanilla < R11.bat > nohup11.out &
#nohup R --slave --vanilla < R12.bat > nohup12.out &
#nohup R --slave --vanilla < R13.bat > nohup13.out &
#nohup R --slave --vanilla < R14.bat > nohup14.out &
#nohup R --slave --vanilla < R15.bat > nohup15.out &
#nohup R --slave --vanilla < R16.bat > nohup16.out &
#nohup R --slave --vanilla < R17.bat > nohup17.out &

files<- list.files(pattern="methyl_P") 
files <- files[c(7,2,3,6,1,4,5)]
methyl_all <-NULL
for (i in c(1:length(files)))
{
load(files[i])
methyl_all <- cbind(methyl_all,methyl[,2] / rowSums(methyl[,2:3]))
}
methyl_all[is.na(methyl_all)] <- 0
pca <- prcomp(scale(methyl_all))
P <- pchisq(scale(pca$x[,1])^2,1,lower.tail=F)
index <- p.adjust(P,"BH")<0.01
methyl[index,1]

#Pcan mRNA expression by pca
x <- read.csv("GSE59525_RPKM_Pcan.txt.gz",sep="\t")
pca1 <- prcomp(scale(2^x[,7:16]))
P1 <- pchisq(scale(pca1$x[,3])^2,1,lower.tail=F)
index1 <- p.adjust(P1,"BH")<0.01
x[index1,1]

#Pcan by HOSVD
require(rTensor)
Z <- array(NA,c(7,10,dim(x)[1]))
for (i in c(1:dim(x)[1]))
{
cat(i," ")
Z[,,i] <- outer(methyl_all[i,],2^x[i,7:16],"*")
}
Z<-apply(Z,c(1,2),scale)
Z <- aperm(Z,c(2,3,1))
HOSVD <- hosvd(as.tensor(Z))
P3 <- pchisq(rowSums(scale(HOSVD$U[[3]][,9:10])^2),2,lower.tail=F)
index3 <- p.adjust(P3,"BH")<0.01
methyl[index3,1]

#--- genome biology ----
x <- read.csv("13059_2012_3057_MOESM9_ESM.CSV",sep=",")
pca <- prcomp(scale(x[,-c(1:2)]))
P <- pchisq(scale(pca$x[,3])^2,1,lower.tail=F)
index <- p.adjust(P,"BH")<0.01
x[index,1]
