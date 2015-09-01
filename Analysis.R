
###############################################################################################
# Preparations
###############################################################################################

# Load packages
library(geomorph)
library(ape)
library(geiger)
library(plyr)

# Read data
data = read.csv("CranialData.csv", header=TRUE, stringsAsFactors=FALSE)
tree = read.nexus("PrimateTree.nex")

# Subset data and phylogeny into males and females
data.m = data[which(data$sex=="M"),]		
data.f = data[which(data$sex=="F"),]	
tree.m = drop.tip(tree, setdiff(tree$tip.label, unique(data.m$genus_species)))
tree.f = drop.tip(tree, setdiff(tree$tip.label, unique(data.f$genus_species)))

# Prepare data for geomorph
array.m = c()
array.f = c()
for (i in unique(data.m$genus_species)) {
  x = data.m[data.m$genus_species==i,]
  array.m = c(array.m, x$x, x$y, x$z)
  }
for (i in unique(data.f$genus_species)) {
  x = data.f[data.f$genus_species==i,] 
  array.f = c(array.f, x$x, x$y, x$z)
  }
landmarks = c("Rhinion", "Nasion", "Bregma", "Medial orbit border", "Lateral orbit border", "Orbitale superior", "Orbitale inferior", "Pterion", "Zygion", "Ectomolare", "Euryon", "Prosthion", "Lambda", "Inion", "Opisthion", "Basion", "Sphenobasion", "Alveolon")
array.m = array(array.m, dim=c(18,3,64), dimnames=list(landmarks, c("x", "y", "z"), unique(data.m$genus_species)))
array.f = array(array.f, dim=c(18,3,61), dimnames=list(landmarks, c("x", "y", "z"), unique(data.f$genus_species)))

###############################################################################################
# 
###############################################################################################