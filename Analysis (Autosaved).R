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
library(adephylo)

# Load functions
source("ProcD.pgls2.R")
source("plotWireframe.R")

# Define links for wireframe diagrams
a = c(1,1,1,2,2,3,13,13,15,16,17,18,4,6,5,4,6,8,8,14,3,9,9,8,10,3,6,6,5)
b = c(2,12,7,4,3,13,14,15,16,17,18,12,6,5,7,7,8,3,11,11,11,11,7,10,12,2,3,7,4)
skull = matrix(c(a,b), length(a), 2)

# Set parameters for wireframe diagrams
#pp <- par3d(no.readonly=TRUE)
#dput(pp, file="figures/wireframe_params/modulesParams.R", control = "all")
pam = dget("figures/wireframe_params/allometryParamsM.R")
paf = dget("figures/wireframe_params/allometryParamsF.R")
pdm = dget("figures/wireframe_params/dimorphismParamsM.R")
pdf = dget("figures/wireframe_params/dimorphismParamsF.R")
pdm2 = dget("figures/wireframe_params/dimorphism2ParamsM.R")
pdf2 = dget("figures/wireframe_params/dimorphism2ParamsF.R")
pm = dget("figures/wireframe_params/modulesParams.R")

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
colorcoding = ddply(landmarks, .(suborder, superfamily, genus_species), function(x) {
  subcol = 0
  supercol = 0
  if (x[1,"suborder"] == "Strepsirrhini") {subcol = "purple"}
  if (x[1,"suborder"] == "Haplorhini") {subcol = "darkgreen"}
  if (x[1,"superfamily"] == "Lorisoidea") {supercol = "black"}
  if (x[1,"superfamily"] == "Lemuroidea") {supercol = "purple"}
  if (x[1,"superfamily"] == "Tarsioidea") {supercol = "blue"}
  if (x[1,"superfamily"] == "Ceboidea") {supercol = "darkgreen"}
  if (x[1,"superfamily"] == "Cercopithecoidea") {supercol = "orange"}
  if (x[1,"superfamily"] == "Hominoidea") {supercol = "red"}
  return(data.frame(subcol=subcol, supercol=supercol))
})
subcol.m = as.character(colorcoding$subcol[match(tree.m$tip.label, colorcoding$genus_species)])
supercol.m = as.character(colorcoding$supercol[match(tree.m$tip.label, colorcoding$genus_species)])
layout(matrix(1:2,1,2))
par(mar=c(1,1,1,2))
plot(tree.m, tip.color=subcol.m, cex=0.35, label.offset=0.5)
mtext("Suborders")
legend("bottomleft", legend=c("Strepsirrhini", "Haplorhini"), fill=c("purple", "darkgreen"), cex=0.4)
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
# males
allometry.m2 = procD.pgls2(coords.m ~ log(csize.m), tree.m)
csize.model.m = coefficients(allometry.m2$model)
min.csize.m = log(csize.m[which.min(csize.m)])
max.csize.m = log(csize.m[which.max(csize.m)])
pred.max.csize.m <- pred.min.csize.m <- c()
for (i in 1:ncol(csize.model.m)) {
  pred.max.csize.m[i] = csize.model.m[1,i] + max.csize.m*csize.model.m[2,i]
  pred.min.csize.m[i] = csize.model.m[1,i] + min.csize.m*csize.model.m[2,i]
}
array.max.csize.m = array(pred.max.csize.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.csize.m = array(pred.min.csize.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.csize.m, skull, params=pam)
plotWireframe(array.min.csize.m, skull, params=pam)
# females
allometry.f2 = procD.pgls2(coords.f ~ log(csize.f), tree.f)
csize.model.f = coefficients(allometry.f2$model)
min.csize.f = log(csize.f[which.min(csize.f)])
max.csize.f = log(csize.f[which.max(csize.f)])
pred.max.csize.f <- pred.min.csize.f <- c()
for (i in 1:ncol(csize.model.f)) {
  pred.max.csize.f[i] = csize.model.f[1,i] + max.csize.f*csize.model.f[2,i]
  pred.min.csize.f[i] = csize.model.f[1,i] + min.csize.f*csize.model.f[2,i]
}
array.max.csize.f = array(pred.max.csize.f, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.csize.f = array(pred.min.csize.f, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.csize.f, skull, params=paf)
plotWireframe(array.min.csize.f, skull, params=paf)

# Get residuals from allometry analysis
resid.coords = function(coords, csize, coefs) {
	dm = dim(coords)
	dn = dimnames(coords)
	for (i in 1:length(csize)) {
		pred = c()
		for (ii in 1:ncol(coefs)) {pred[ii] = coefs[1,ii] + csize[i]*coefs[2,ii]}
		res = two.d.array(coords)[i,] - pred
		coords[,,i] = matrix(res, dm[1:2], byrow=TRUE, dimnames=dn[1:2])
	}
	return(coords)	
}
resid.m = resid.coords(coords.m, log(csize.m), coefficients(allometry.m2$model))
resid.f = resid.coords(coords.f, log(csize.f), coefficients(allometry.f2$model))

# Phylomorphospace plots for residuals of allometry analysis
quartz()
layout(matrix(1:2,1,2))
plotGMPhyloMorphoSpace(tree.m, resid.m, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.3))
mtext("Males", line=1)
plotGMPhyloMorphoSpace(tree.f, resid.f, plot.param=list(t.cex=0.3, n.cex=0.3, lwd=0.3, txt.cex=0.3))
mtext("Females", line=1)

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

# Wireframes of dimorphism for males and females (controlling for centroid size)
# males
sexselect.m2 = procD.pgls2(coords.m ~ dimorphism.m + log(csize.m), tree.m)
dimorph.model.m = coefficients(sexselect.m2$model)
min.dimorph.m = dimorphism.m[which.min(dimorphism.m)]
max.dimorph.m = dimorphism.m[which.max(dimorphism.m)]
pred.max.dimorph.m <- pred.min.dimorph.m <- c()
for (i in 1:ncol(dimorph.model.m)) {
  pred.max.dimorph.m[i] = dimorph.model.m[1,i] + max.dimorph.m*dimorph.model.m[2,i] + mean(log(csize.m))*dimorph.model.m[3,i]
  pred.min.dimorph.m[i] = dimorph.model.m[1,i] + min.dimorph.m*dimorph.model.m[2,i] + mean(log(csize.m))*dimorph.model.m[3,i]
}
array.max.dimorph.m = array(pred.max.dimorph.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.dimorph.m = array(pred.min.dimorph.m, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.dimorph.m, skull, params=pdm)
plotWireframe(array.min.dimorph.m, skull, params=pdm)
# females
sexselect.f2 = procD.pgls2(coords.f ~ dimorphism.f + log(csize.f), tree.f)
dimorph.model.f = coefficients(sexselect.f2$model)
min.dimorph.f = dimorphism.f[which.min(dimorphism.f)]
max.dimorph.f = dimorphism.f[which.max(dimorphism.f)]
pred.max.dimorph.f <- pred.min.dimorph.f <- c()
for (i in 1:ncol(dimorph.model.f)) {
  pred.max.dimorph.f[i] = dimorph.model.f[1,i] + max.dimorph.f*dimorph.model.f[2,i] + mean(log(csize.f))*dimorph.model.f[3,i]
  pred.min.dimorph.f[i] = dimorph.model.f[1,i] + min.dimorph.f*dimorph.model.f[2,i] + mean(log(csize.f))*dimorph.model.f[3,i]
}
array.max.dimorph.f = array(pred.max.dimorph.f, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.dimorph.f = array(pred.min.dimorph.f, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.dimorph.f, skull, params=pdf)
plotWireframe(array.min.dimorph.f, skull, params=pdf)

# Wireframes of allometry for males and females (controlling for dimorphism)
# males
pred.max.dimorph.m2 <- pred.min.dimorph.m2 <- c()
for (i in 1:ncol(dimorph.model.m)) {
  pred.max.dimorph.m2[i] = dimorph.model.m[1,i] + mean(dimorphism.m)*dimorph.model.m[2,i] + max.csize.m*dimorph.model.m[3,i]
  pred.min.dimorph.m2[i] = dimorph.model.m[1,i] + mean(dimorphism.m)*dimorph.model.m[2,i] + min.csize.m*dimorph.model.m[3,i]
}
array.max.dimorph.m2 = array(pred.max.dimorph.m2, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.dimorph.m2 = array(pred.min.dimorph.m2, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.dimorph.m2, skull, params=pdm2)
plotWireframe(array.min.dimorph.m2, skull, params=pdm2)
# females
pred.max.dimorph.f2 <- pred.min.dimorph.f2 <- c()
for (i in 1:ncol(dimorph.model.f)) {
  pred.max.dimorph.f2[i] = dimorph.model.f[1,i] + mean(dimorphism.f)*dimorph.model.f[2,i] + max.csize.f*dimorph.model.f[3,i]
  pred.min.dimorph.f2[i] = dimorph.model.f[1,i] + mean(dimorphism.f)*dimorph.model.f[2,i] + min.csize.f*dimorph.model.f[3,i]
}
array.max.dimorph.f2 = array(pred.max.dimorph.f2, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
array.min.dimorph.f2 = array(pred.min.dimorph.f2, dim=c(18,3,1), dimnames=list(landmark.names, c("x", "y", "z")))
plotWireframe(array.max.dimorph.f2, skull, params=pdf2)
plotWireframe(array.min.dimorph.f2, skull, params=pdf2)

###############################################################################################
# Modularity, integration, and rates of cranial shape evolution
###############################################################################################

# Define neurocranium and facial modules
face.names = c("Rhinion", "Nasion", "Medial orbit border", "Lateral orbit border", "Orbitale superior", "Orbitale inferior", "Zygion", "Ectomolare", "Prosthion", "Alveolon")
neuro.names = c("Bregma", "Pterion", "Euryon", "Lambda", "Inion", "Opisthion", "Basion", "Sphenobasion")
face.m = coords.m[face.names,,]
neuro.m = coords.m[neuro.names,,]
face.f = coords.f[face.names,,]
neuro.f = coords.f[neuro.names,,]
partition = c("F", "F", "N", "F", "F", "F", "F", "N", "F", "F", "N", "F", "N", "N", "N", "N", "N", "F")

# Define taxonomic groups to compare
suborder.m = dimnames(A)[[3]]
suborder.m[c(1,2,44:63)] = "Prosimian"
suborder.m[3:43] = "Anthropoid"
cats = dimnames(A)[[3]]
cats[1:length(groups2)] = "NonCat"
cats[4:30] = "Cat"

# Plot wireframe with color coded modules
open3d()
par3d(pm)
rgl.bg(sphere=TRUE, color=c("white"), lit=FALSE, back="fill")
A = coords.m[,,"Cebus_apella"]
Af = coords.m[face.names,,"Cebus_apella"]
An = coords.m[neuro.names,,"Cebus_apella"]
A3d = matrix(A, dim(A)[[1]], dim(A)[[2]], dimnames=dimnames(A)[1:2])
plot3d(A3d, type = "n", col = "red", xlab = "", ylab = "", zlab = "", size = 1, aspect = FALSE, box=FALSE, axes=FALSE)
points3d(Af, color = "red", size = 5)
points3d(An, color = "blue", size = 5)
for (i in 1:nrow(skull)) {segments3d(rbind(A[skull[i,1],], A[skull[i,2],]), lwd = 2, col="black")}

# Define function to compute independent contrasts
pic.shape = function (data, tree) {
  num_landmarks = length(dimnames(data)[[1]])
  num_dimensions = length(dimnames(data)[[2]])
  num_species = length(dimnames(data)[[3]])
  new_array = array(rep(0,length(data)), dim=c(num_landmarks, num_dimensions, num_species-1))
  for (i in 1:num_landmarks) {for (ii in 1:num_dimensions) {new_array[i,ii,] = pic(data[i,ii,], tree)}}
  dimnames(new_array) = list(dimnames(data)[[1]], dimnames(data)[[2]], names(pic(data[1,1,], tree)))
  return(new_array)
}

# Test for modularity (raw coords and size controlled)
pic.coords.m = pic.shape(coords.m, tree.m)
pic.coords.f = pic.shape(coords.f, tree.f)
compare.modular.partitions(pic.coords.m, partition)# p = 0.005, rv=0.58
compare.modular.partitions(pic.coords.f, partition)# p=0.002, rv=0.47
pic.resid.m = pic.shape(resid.m, tree.m)
pic.resid.f = pic.shape(resid.f, tree.f)
compare.modular.partitions(pic.resid.m, partition)# p = 0.002, rv=0.52
compare.modular.partitions(pic.resid.f, partition)# p=0.001, rv=0.45

# Test for integration (raw coords and size controlled)
morphol.integr(pic.coords.m[face.names,,], pic.coords.m[neuro.names,,], method="RV")# p=0.001
morphol.integr(pic.coords.f[face.names,,], pic.coords.f[neuro.names,,], method="RV")# p=0.001
morphol.integr(pic.resid.m[face.names,,], pic.resid.m[neuro.names,,], method="RV")# p=0.001
morphol.integr(pic.resid.f[face.names,,], pic.resid.f[neuro.names,,], method="RV")# p=0.001

# Test for differences in integration/modularity between strepsirrhines and haplorhines (raw coords)
streps.m = coords.m[,,unique(landmarks.m$genus_species[which(landmarks.m$suborder == "Strepsirrhini")])]
streps.m = gpagen(streps.m, ShowPlot = FALSE)
streps.m = array(streps.m$coords, dim=dim(streps.m$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(streps.m$coords)[[3]]))
streps.f = coords.f[,,unique(landmarks.f$genus_species[which(landmarks.f$suborder == "Strepsirrhini")])]
streps.f = gpagen(streps.f, ShowPlot = FALSE)
streps.f = array(streps.f$coords, dim=dim(streps.f$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(streps.f$coords)[[3]]))
haps.m = coords.m[,,unique(landmarks.m$genus_species[which(landmarks.m$suborder == "Haplorhini")])]
haps.m = gpagen(haps.m, ShowPlot = FALSE)
haps.m = array(haps.m$coords, dim=dim(haps.m$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(haps.m$coords)[[3]]))
haps.f = coords.f[,,unique(landmarks.f$genus_species[which(landmarks.f$suborder == "Haplorhini")])]
haps.f = gpagen(haps.f, ShowPlot = FALSE)
haps.f = array(haps.f$coords, dim=dim(haps.f$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(haps.f$coords)[[3]]))
tree.strep.m = drop.tip(tree.m, setdiff(tree.m$tip.label, dimnames(streps.m)[[3]]))
tree.strep.f = drop.tip(tree.f, setdiff(tree.f$tip.label, dimnames(streps.f)[[3]]))
tree.haps.m = drop.tip(tree.m, setdiff(tree.m$tip.label, dimnames(haps.m)[[3]]))
tree.haps.f = drop.tip(tree.f, setdiff(tree.f$tip.label, dimnames(haps.f)[[3]]))
pic.streps.m = pic.shape(streps.m, tree.strep.m)
pic.streps.f = pic.shape(streps.f, tree.strep.f)
pic.haps.m = pic.shape(haps.m, tree.haps.m)
pic.haps.f = pic.shape(haps.f, tree.haps.f)
compare.modular.partitions(pic.streps.m, partition)#p=0.5, rv=0.66
compare.modular.partitions(pic.haps.m, partition)#p=0.001, rv=0.6
compare.modular.partitions(pic.streps.f, partition)#p=0.13, rv=0.65
compare.modular.partitions(pic.haps.f, partition)#p=0.001, rv=0.5
morphol.integr(pic.streps.m[neuro.names,,], pic.streps.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.haps.m[neuro.names,,], pic.haps.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.streps.f[neuro.names,,], pic.streps.f[face.names,,], method="RV")#p=0.002
morphol.integr(pic.haps.f[neuro.names,,], pic.haps.f[face.names,,], method="RV")#p=0.001

# Test for differences in integration/modularity between strepsirrhines and haplorhines (size controlled)
streps.m2 = resid.m[,,unique(landmarks.m$genus_species[which(landmarks.m$suborder == "Strepsirrhini")])]
streps.m2 = gpagen(streps.m2, ShowPlot = FALSE)
streps.m2 = array(streps.m2$coords, dim=dim(streps.m2$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(streps.m2$coords)[[3]]))
streps.f2 = resid.f[,,unique(landmarks.f$genus_species[which(landmarks.f$suborder == "Strepsirrhini")])]
streps.f2 = gpagen(streps.f2, ShowPlot = FALSE)
streps.f2 = array(streps.f2$coords, dim=dim(streps.f2$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(streps.f2$coords)[[3]]))
haps.m2 = resid.m[,,unique(landmarks.m$genus_species[which(landmarks.m$suborder == "Haplorhini")])]
haps.m2 = gpagen(haps.m2, ShowPlot = FALSE)
haps.m2 = array(haps.m2$coords, dim=dim(haps.m2$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(haps.m2$coords)[[3]]))
haps.f2 = resid.f[,,unique(landmarks.f$genus_species[which(landmarks.f$suborder == "Haplorhini")])]
haps.f2 = gpagen(haps.f2, ShowPlot = FALSE)
haps.f2 = array(haps.f2$coords, dim=dim(haps.f2$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(haps.f2$coords)[[3]]))
pic.streps.m2 = pic.shape(streps.m2, tree.strep.m)
pic.streps.f2 = pic.shape(streps.f2, tree.strep.f)
pic.haps.m2 = pic.shape(haps.m2, tree.haps.m)
pic.haps.f2 = pic.shape(haps.f2, tree.haps.f)
compare.modular.partitions(pic.streps.m2, partition)#p=0.45, rv=0.62
compare.modular.partitions(pic.haps.m2, partition)#p=0.001, rv=0.47
compare.modular.partitions(pic.streps.f2, partition)#p=0.24, rv=0.65
compare.modular.partitions(pic.haps.f2, partition)#p=0.002, rv=0.51
# significant difference between groups?
compare.modular.groups(pic.resid.m, partition, suborder.m, ShowPlot=TRUE)
morphol.integr(pic.streps.m[neuro.names,,], pic.streps.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.haps.m[neuro.names,,], pic.haps.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.streps.f[neuro.names,,], pic.streps.f[face.names,,], method="RV")#p=0.002
morphol.integr(pic.haps.f[neuro.names,,], pic.haps.f[face.names,,], method="RV")#p=0.001

# Test for different rates of evolution in strepsirrhines and haplorhines
factors.m = factor(unique(landmarks.m[,c("genus_species", "suborder")])[,2])
factors.f = factor(unique(landmarks.f[,c("genus_species", "suborder")])[,2])
names(factors.m) = unique(landmarks.m$genus_species)
names(factors.f) = unique(landmarks.f$genus_species)
rates.m = compare.evol.rates(tree.m, coords.m, factors.m, ShowPlot=FALSE)#p=0.001
2.286179e-05/8.129630e-06 
rates.f = compare.evol.rates(tree.f, coords.f, factors.f, ShowPlot=FALSE)#p=0.001
1.645051e-05/8.085462e-06 

# Test for differences in integration/modularity between catarrhini and platyrrhini
cats.m = coords.m[,,unique(landmarks.m$genus_species[which(landmarks.m$superfamily %in% c("Cercopithecoidea", "Hominoidea"))])]
cats.m = gpagen(cats.m, ShowPlot = FALSE)
cats.m = array(cats.m$coords, dim=dim(cats.m$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(cats.m$coords)[[3]]))
cats.f = coords.f[,,unique(landmarks.f$genus_species[which(landmarks.f$superfamily %in% c("Cercopithecoidea", "Hominoidea"))])]
cats.f = gpagen(cats.f, ShowPlot = FALSE)
cats.f = array(cats.f$coords, dim=dim(cats.f$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(cats.f$coords)[[3]]))
plats.m = coords.m[,,unique(landmarks.m$genus_species[which(landmarks.m$superfamily == "Ceboidea")])]
plats.m = gpagen(plats.m, ShowPlot = FALSE)
plats.m = array(plats.m$coords, dim=dim(plats.m$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(plats.m$coords)[[3]]))
plats.f = coords.f[,,unique(landmarks.f$genus_species[which(landmarks.f$superfamily == "Ceboidea")])]
plats.f = gpagen(plats.f, ShowPlot = FALSE)
plats.f = array(plats.f$coords, dim=dim(plats.f$coords), dimnames=list(landmark.names, c("x", "y", "z"), dimnames(plats.f$coords)[[3]]))
tree.cats.m = drop.tip(tree.m, setdiff(tree.m$tip.label, dimnames(cats.m)[[3]]))
tree.cats.f = drop.tip(tree.f, setdiff(tree.f$tip.label, dimnames(cats.f)[[3]]))
tree.plats.m = drop.tip(tree.m, setdiff(tree.m$tip.label, dimnames(plats.m)[[3]]))
tree.plats.f = drop.tip(tree.f, setdiff(tree.f$tip.label, dimnames(plats.f)[[3]]))
pic.cats.m = pic.shape(cats.m, tree.cats.m)
pic.cats.f = pic.shape(cats.f, tree.cats.f)
pic.plats.m = pic.shape(plats.m, tree.plats.m)
pic.plats.f = pic.shape(plats.f, tree.plats.f)
compare.modular.partitions(pic.cats.m, partition)#p=0.001, rv=0.64
compare.modular.partitions(pic.plats.m, partition)#p=0.586, rv=0.77
compare.modular.partitions(pic.cats.f, partition)#p=0.001, rv=0.53
compare.modular.partitions(pic.plats.f, partition)#p=0.012, rv=0.54
morphol.integr(pic.cats.m[neuro.names,,], pic.cats.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.plats.m[neuro.names,,], pic.plats.m[face.names,,], method="RV")#p=0.001
morphol.integr(pic.cats.f[neuro.names,,], pic.cats.f[face.names,,], method="RV")#p=0.001
morphol.integr(pic.plats.f[neuro.names,,], pic.plats.f[face.names,,], method="RV")#p=0.017

# Test for different rates of evolution in catarrhini and platyrrhini
tree.anth.m = drop.tip(tree.m, setdiff(tree.m$tip.label, c(dimnames(cats.m)[[3]], dimnames(plats.m)[[3]])))
tree.anth.f = drop.tip(tree.f, setdiff(tree.f$tip.label, c(dimnames(cats.f)[[3]], dimnames(plats.f)[[3]])))
coords.anth.m = coords.m[,,unique(landmarks.m$genus_species[which(landmarks.m$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea"))])]
coords.anth.f = coords.f[,,unique(landmarks.f$genus_species[which(landmarks.f$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea"))])]
factors.anth.m = unique(landmarks.m[which(landmarks.m$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea")),c("genus_species", "superfamily")])[,2]
factors.anth.m[which(factors.anth.m  %in% c("Cercopithecoidea", "Hominoidea"))] = "Catarrhini"
factors.anth.m[which(factors.anth.m  == "Ceboidea")] = "Platyrrhini"
factors.anth.m = factor(factors.anth.m)
names(factors.anth.m) = unique(landmarks.m[which(landmarks.m$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea")),c("genus_species", "superfamily")])[,1]
factors.anth.f = unique(landmarks.f[which(landmarks.f$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea")),c("genus_species", "superfamily")])[,2]
factors.anth.f[which(factors.anth.f  %in% c("Cercopithecoidea", "Hominoidea"))] = "Catarrhini"
factors.anth.f[which(factors.anth.f  == "Ceboidea")] = "Platyrrhini"
factors.anth.f = factor(factors.anth.f)
names(factors.anth.f) = unique(landmarks.f[which(landmarks.f$superfamily %in% c("Cercopithecoidea", "Hominoidea", "Ceboidea")),c("genus_species", "superfamily")])[,1]
rates.anth.m = compare.evol.rates(tree.anth.m, coords.anth.m, factors.anth.m, ShowPlot=FALSE)#p=0.026
2.65970e-05/1.63861e-05 
rates.anth.f = compare.evol.rates(tree.anth.f, coords.anth.f, factors.anth.f, ShowPlot=FALSE)#p=0.015
1.898842e-05/1.234976e-05









