# Set working directory
setwd("~/Desktop/GitHub/Primate_cranial_evolution")

# In case user is running windows and can't use quartz 
if (.Platform$OS.type=="windows") {quartz<-function() windows()}

# Load packages
library(geomorph)
library(ape)
library(geiger)
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












