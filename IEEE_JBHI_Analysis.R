library(bio3d)
library(ggplot2)
library(ggpubr)
# Read multi-model PDB file
WT <- read.pdb(file.choose(), multi=TRUE)
MT <- read.pdb(file.choose(), multi=TRUE)
MRT <- read.pdb(file.choose(),multi = TRUE)
C797S<- read.pdb(file.choose(),multi = TRUE)

## RMSD Analysis ##
xyz <- pdbfit(WT)
ca.inds <- atom.select(WT, elety="CA")
rd <- rmsd(WT$xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
rdw <- rd
plot(rdw, main="RMSD of EGFR and its mutants", ylab="RMSD in A", xlab="Frame Number", typ="l", xlim =(c(1,1000)), ylim = (c(0, 5)))

xyz1 <- pdbfit(MT)
ca.inds1 <- atom.select(MT, elety="CA")
rd1 <- rmsd(MT$xyz[1,ca.inds1$xyz], xyz[,ca.inds1$xyz])
plot(rd1, main="RMSD of EGFR and its mutants", ylab="RMSD in A", xlab="Frame Number", typ="l", xlim =(c(1,1000)), ylim = (c(0, 6)))

xyz2 <- pdbfit(MRT)
ca.inds2 <- atom.select(MRT, elety="CA")
rd2 <- rmsd(MRT$xyz[1,ca.inds2$xyz], xyz[,ca.inds2$xyz])
plot(rd2, main="RMSD of EGFR and its mutants", ylab="RMSD in A", xlab="Frame Number", typ="l", xlim =(c(1,1000)), ylim = (c(0, 6)))

xyz3 <- pdbfit(C797S)
ca.inds3 <- atom.select(C797S, elety="CA")
rd3 <- rmsd(C797S$xyz[1,ca.inds3$xyz], xyz[,ca.inds3$xyz])
rdw3 <- rd
plot(rdw3, main="RMSD of EGFR and its mutants", ylab="RMSD in A", xlab="Frame Number", typ="l", xlim =(c(1,1000)), ylim = (c(0, 5)))



##Plots of multiple rmsds on a single figure %%%
plot(rd, main="RMSD of EGFR and its mutants", ylab="RMSD in A", xlab="Frame Number", typ="l", xlim =(c(1,1000)), ylim = (c(0, 10)))
lines(rd1, typ ="l",col="blue")
lines(rd2,,col="red")
lines(rd3,,col="green")
lines(rd4,,col="orange")

grid(10,10)

legend("bottomright",
       c("WT","L858R","LT790M","C797S"),
       fill=c("black","blue", "red","green"))


# Principal Component Analysis (PCA)
pc <- pca.xyz(xyz[,ca.inds3$xyz])
# Plot PCA
plot(pc, col=bwr.colors(nrow(xyz)))

## Contact map ##
ca.inds <- atom.select(WT, "calpha")
a <- WT$xyz[1,ca.inds$xyz]
ref.cont <- cmap(a, dcut=10, scut=1)
count3 <- length(which(ref.cont == 1))
par(mfcol = c(1,1))
plot.cmap(ref.cont, ylab="Residue Position", xlab="Residue Position", main ="Contact map Mutated and resistive type EGFR Kinase Domain", type = "l", col = "red")
count3
ref.cont
sum.cont <- NULL 
for(i in 1:nrow(dw7_avg)) {
  ## Contact map for frame 'i' 
  cont <- cmap(dw7[i,inds$xyz], dcut=6, scut=3)
  ## Product with reference 
  prod.cont <- ref.cont.cont 
  sum.cont <- c(sum.cont, sum(prod.cont,na.rm=TRUE))
}
plot(sum.cont, typ="l")
