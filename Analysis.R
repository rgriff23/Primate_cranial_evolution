# Set working directory
setwd("~/Desktop/GitHub/Primate_cranial_evolution")

# In case user is running windows and can't use quartz 
if (.Platform$OS.type=="windows") {quartz<-function() windows()}

# Load packages
library(geomorph)
library(ape)
library(geiger)
library(caper)
library(plyr)

# Set seed since some significance tests are based on permutations
set.seed(922015)

###############################################################################################
# Preparations
###############################################################################################

# Read data
landmarks = read.csv("./data/CranialData.csv", header=TRUE, stringsAsFactors=FALSE)
tree = read.nexus("./data/PrimateTree.nex")
ecology2 = read.csv("./data/EcologyData_2state.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
ecology3 = read.csv("./data/EcologyData_3state.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)
bodymass = read.csv("./data/BodyMass.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

# Subset landmarks and phylogeny into males and females
landmarks.m = landmarks[which(landmarks$sex=="M"),]		
landmarks.f = landmarks[which(landmarks$sex=="F"),]	
tree.m = drop.tip(tree, setdiff(tree$tip.label, unique(landmarks.m$genus_species)))
tree.f = drop.tip(tree, setdiff(tree$tip.label, unique(landmarks.f$genus_species)))

# Format landmark data as 3D arrays for geomorph
array.m = c()
array.f = c()
for (i in unique(landmarks.m$genus_species)) {
  x = landmarks.m[landmarks.m$genus_species==i,]
  array.m = c(array.m, x$x, x$y, x$z)
  }
for (i in unique(landmarks.f$genus_species)) {
  x = landmarks.f[landmarks.f$genus_species==i,] 
  array.f = c(array.f, x$x, x$y, x$z)
  }
landmark.names = c("Rhinion", "Nasion", "Bregma", "Medial orbit border", "Lateral orbit border", "Orbitale superior", "Orbitale inferior", "Pterion", "Zygion", "Ectomolare", "Euryon", "Prosthion", "Lambda", "Inion", "Opisthion", "Basion", "Sphenobasion", "Alveolon")
array.m = array(array.m, dim=c(18,3,64), dimnames=list(landmark.names, c("x", "y", "z"), unique(landmarks.m$genus_species)))
array.f = array(array.f, dim=c(18,3,61), dimnames=list(landmark.names, c("x", "y", "z"), unique(landmarks.f$genus_species)))

# View species and phylogeny
quartz()
colorcoding = ddply(landmarks, .(infraorder, superfamily, genus_species), function(x) {
  infracol = 0
  supercol = 0
  if (x[1,"infraorder"] == "Lemuriformes") {infracol = "purple"}
  if (x[1,"infraorder"] == "Tarsiiformes") {infracol = "blue"}
  if (x[1,"infraorder"] == "Simiiformes") {infracol = "darkgreen"}
  if (x[1,"superfamily"] == "Lorisoidea") {supercol = "black"}
  if (x[1,"superfamily"] == "Lemuroidea") {supercol = "purple"}
  if (x[1,"superfamily"] == "Tarsioidea") {supercol = "blue"}
  if (x[1,"superfamily"] == "Ceboidea") {supercol = "darkgreen"}
  if (x[1,"superfamily"] == "Cercopithecoidea") {supercol = "orange"}
  if (x[1,"superfamily"] == "Hominoidea") {supercol = "red"}
  return(data.frame(infracol=infracol, supercol=supercol))
})
infracol.m = as.character(colorcoding$infracol[match(tree.m$tip.label, colorcoding$genus_species)])
supercol.m = as.character(colorcoding$supercol[match(tree.m$tip.label, colorcoding$genus_species)])
layout(matrix(1:2,1,2))
par(mar=c(1,1,1,2))
plot(tree.m, tip.color=infracol.m, cex=0.35, label.offset=0.5)
mtext("Infraorders")
legend("bottomleft", legend=c("Lemuriformes", "Tarsiiformes", "Simiiformes"), fill=c("purple", "blue", "darkgreen"), cex=0.4)
plot(tree.m, tip.color=supercol.m, cex=0.35, label.offset=0.5)
mtext("Superfamilies")
legend("bottomleft", legend=c("Lorisoidea", "Lemuroidea", "Tarsioidea", "Ceboidea", "Cercopithecoidea", "Hominoidea"), fill=c("black", "purple", "blue", "darkgreen", "orange", "red"), cex=0.4)

###############################################################################################
# Generalized Procrustes Alignment and testing for phylogenetic signal
###############################################################################################

# Generalized Procrustes Alignment (GPA)
gpa.m = gpagen(array.m, ShowPlot=FALSE)
gpa.f = gpagen(array.f, ShowPlot=FALSE)

# Store GPA coordinates and centroid sizes separately
coords.m = array(gpa.m$coords, dim=c(18,3,64), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(gpa.m$coords)[[3]]))
coords.f = array(gpa.f$coords, dim=c(18,3,61), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(gpa.f$coords)[[3]]))
csize.m = gpa.m$Csize
csize.f = gpa.f$Csize

# Testing for phylogenetic signal with multivariate Blomberg's K
phy.kmult.m = physignal(tree.m, coords.m, method="Kmult", ShowPlot=FALSE)
phy.kmult.f = physignal(tree.f, coords.f, method="Kmult", ShowPlot=FALSE)

# Testing for phylogenetic signal with Klingenberg and Marugan-Lobon's SCP method
phy.ssc.m = physignal(tree.m, coords.m, method="SSC", ShowPlot=FALSE)
phy.ssc.f = physignal(tree.f, coords.f, method="SSC", ShowPlot=FALSE)

# Phylomorphospace figures
quartz()
layout(matrix(1:2,1,2))
plotGMPhyloMorphoSpace(tree.m, coords.m, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.3))
mtext("Males", line=1)
plotGMPhyloMorphoSpace(tree.f, coords.f, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.3))
mtext("Females", line=1)

###############################################################################################
# Allometry of cranial shape
###############################################################################################

# Run regressions of shape vs size
allometry.m = procD.pgls(coords.m ~ log(csize.m), tree.m)
allometry.f = procD.pgls(coords.f ~ log(csize.f), tree.f)

# Wireframes for allometry for males and females
source("ProcD.pgls2.R")
allometry.m2 = procD.pgls2(coords.m ~ log(csize.m), tree.m)
x = data.frame(c(log(csize.m["Microcebus_murinus"]), log(csize.m["Gorilla_gorilla"])))
names(x) = "x.newlog(csize.m)"
test = predict(allometry.m2$model)

csize.model.m = coefficients(allometry.m2$model)
min.csize.m = log(csize.m[which.min(csize.m)])
max.csize.m = log(csize.m[which.max(csize.m)])
pred.max.csize.m <- pred.min.csize.m <- c()
for (i in 1:ncol(csize.model.m)) {
  pred.max.csize.m[i] = csize.model.m[1,i] + max.csize.m*csize.model.m[2,i]
  pred.min.csize.m[i] = csize.model.m[1,i] + min.csize.m*csize.model.m[2,i]
}
a = c(1,1,1,2,2,3,13,13,15,16,17,18,4,6,5,4,6,8,8,14,3,9,9,8,10)
b = c(2,12,7,4,3,13,14,15,16,17,18,12,6,5,7,7,8,3,11,11,11,11,7,10,12)
skull = matrix(c(a,b), 25, 2)
array.max.csize.m = array(pred.max.csize.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")), links=skull)
array.min.csize.m = array(pred.min.csize.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")), links=skull)
plotAllSpecimens(array.max.csize.m)
plotAllSpecimens(array.min.csize.m)

###############################################################################################
# Ecomorphology of cranial shape
###############################################################################################

# Prepare male ecology data
ecology2.m = ecology2[names(csize.m),]
activity2.m = setNames(ecology2.m$activity, rownames(ecology2.m))
diet2.m = setNames(ecology2.m$diet, rownames(ecology2.m))
locomotion2.m = setNames(ecology2.m$locomotion, rownames(ecology2.m))

# Prepare female ecology data
ecology2.f = ecology2[names(csize.f),]
activity2.f = setNames(ecology2.f$activity, rownames(ecology2.f))
diet2.f = setNames(ecology2.f$diet, rownames(ecology2.f))
locomotion2.f = setNames(ecology2.f$locomotion, rownames(ecology2.f))

# Plot ecology data
quartz()
layout(matrix(1:3,1,3))
par(mar=c(1,1,1,2))
act_col = c("orange", "black")
diet_col = c("hotpink", "purple", "blue")
loc_col = c("saddlebrown", "green4")
plot(tree, tip.color=act_col[ecology2$activity][match(tree$tip.label, rownames(ecology2))], cex=0.7, label.offset=0.5)
legend("bottomleft", legend=c("Diurnal/Cathemeral", "Nocturnal"), fill=act_col, cex=0.7)
plot(tree, tip.color=diet_col[ecology2$diet][match(tree$tip.label, rownames(ecology2))], cex=0.7, label.offset=0.5)
legend("bottomleft", legend=c("Frugivore", "Omnivore", "Insectivore/Folivore"), fill=diet_col, cex=0.7)
plot(tree, tip.color=loc_col[ecology2$locomotion][match(tree$tip.label, rownames(ecology2))], cex=0.7, label.offset=0.5)
legend("bottomleft", legend=c("Pronograde", "Orthograde"), fill=loc_col, cex=0.7)
par(mar=c(5,4,4,2)+0.1)

# Run regressions of shape vs ecology
ecomorph.m = procD.pgls(coords.m ~ activity2.m + diet2.m + locomotion2.m + log(csize.m), tree.m, RRPP=TRUE)
ecomorph.f = procD.pgls(coords.f ~ activity2.f + diet2.f + locomotion2.f + log(csize.f), tree.f, RRPP=TRUE)

# Run regressions of shape vs ecology with Pagel's lambda = 0
ecomorph.m0 = procD.pgls(coords.m ~ activity2.m + diet2.m + locomotion2.m + log(csize.m), rescale(tree.m, "lambda", 0), RRPP=TRUE)
ecomorph.f0 = procD.pgls(coords.f ~ activity2.f + diet2.f + locomotion2.f + log(csize.f), rescale(tree.f, "lambda", 0), RRPP=TRUE)

# Wireframes for nocturnality

###############################################################################################
# Sexual dimorphism and cranial shape
###############################################################################################

# Prepare sexual dimorphism data
dimorphism.m = bodymass[names(csize.m),]
dimorphism.m = setNames(log(dimorphism.m$male_bodymass) - log(dimorphism.m$female_bodymass), names(csize.m))
dimorphism.f = bodymass[names(csize.f),]
dimorphism.f = setNames(log(dimorphism.f$male_bodymass) - log(dimorphism.f$female_bodymass), names(csize.f))

# Variance inflation factors for body mass and sexual dimorphism
compdat.m = comparative.data(tree.m, data.frame(dimorphism.m, csize.m, species=names(csize.m)), species)
compdat.f = comparative.data(tree.f, data.frame(dimorphism.f, csize.f, species=names(csize.f)), species)
lm.bm.sd.m = pgls(dimorphism.m ~ log(csize.m), compdat.m)
lm.bm.sd.f = pgls(dimorphism.f ~ log(csize.f), compdat.f)
vif.bm.sd.m = 1/(1 - summary(lm.bm.sd.m)$adj.r.squared)
vif.bm.sd.f = 1/(1 - summary(lm.bm.sd.f)$adj.r.squared)
sqrt(vif.bm.sd.m)
sqrt(vif.bm.sd.f)

# Run regressions of shape vs. dimorphism
sexselect.m = procD.pgls(coords.m ~ dimorphism.m + log(csize.m), tree.m, RRPP=TRUE)
sexselect.f = procD.pgls(coords.f ~ dimorphism.f + log(csize.f), tree.f, RRPP=TRUE)

# Wireframes of dimorphism for males and females

###############################################################################################
# Sexual dimorphism and cranial shape divergence
###############################################################################################

# Compute dimorphism PIC ancestral states
dimorphism.m.ace = ace(dimorphism.m, tree.m, method="pic")$ace
dimorphism.f.ace = ace(dimorphism.f, tree.f, method="pic")$ace

# Compute centroid size PICs
csize.pic.m = pic(log(csize.m), tree.m)
csize.pic.f = pic(log(csize.f), tree.f)

# Compute sum of shape PICs
mat.m = matrix(coords.m, dim(coords.m)[3], prod(dim(coords.m)[1:2]), byrow=TRUE, dimnames=list(dimnames(coords.m)[[3]]))
mat.f = matrix(coords.f, dim(coords.f)[3], prod(dim(coords.f)[1:2]), byrow=TRUE, dimnames=list(dimnames(coords.f)[[3]]))
shape.pic.m = c()
shape.pic.f = c()
for (i in 1:ncol(mat.m)) {
  shape.pic.m = cbind(shape.pic.m, pic(mat.m[,i], tree.m))
  shape.pic.f = cbind(shape.pic.f, pic(mat.f[,i], tree.f))
  
}
sum.pic.m = rowSums(abs(shape.pic.m))
sum.pic.f = rowSums(abs(shape.pic.f))

# Regress shape disparity vs. dimorphism and centroid size
dimorphism.lm.m = lm(log(sum.pic.m) ~ dimorphism.m.ace + csize.pic.m)
dimorphism.lm.f = lm(log(sum.pic.f) ~ dimorphism.f.ace + csize.pic.f)
summary(dimorphism.lm.m)$coef
summary(dimorphism.lm.f)$coef

# Visualize results with scatterplots
quartz()
layout(matrix(1:2, 1, 2))
pch.m = c(rep(1, 44), rep(3, 19))
pch.f = c(rep(1, 43), rep(3, 17))
plot(log(sum.pic.m) ~ dimorphism.m.ace, ylab="Cranial shape divergence", xlab="Sexual size dimorphism", main="Males", pch=pch.m)
abline(a=summary(dimorphism.lm.m)$coef[1,"Estimate"], b=summary(dimorphism.lm.m)$coef[2,"Estimate"])
par(xpd=TRUE)
legend("topleft", inset=c(-0.25, -0.25), legend=c("Haplorhines", "Strepsirrhines"), pch=c(1,3), cex=0.6)
par(xpd=FALSE)
plot(log(sum.pic.f) ~ dimorphism.f.ace, ylab="", xlab="Sexual size dimorphism", main="Females", pch=pch.f)
abline(a=summary(dimorphism.lm.f)$coef[1,"Estimate"], b=summary(dimorphism.lm.f)$coef[2,"Estimate"])

# Visualize magnitude of cranial shape changes across phylogeny based on PIC
quartz()
layout(matrix(1:2,1,2))
par(mar=c(0,0,2,0))
plot(tree.m, cex=0.5)
nodelabels(pch=22, cex=sum.pic.m*4, col="blue")
mtext("Males")
plot(tree.f, cex=0.5)
nodelabels(pch=22, cex=sum.pic.f*4, col="red", bg="pink")
mtext("Females")

###############################################################################################
# Modularity, integration, and rates of cranial shape evolution
###############################################################################################









