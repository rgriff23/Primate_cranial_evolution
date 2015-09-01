# Primate cranial shape evolution: Combining geometric morphometrics and phylogenetic comparative methods

[Randi H. Griffin]()

* The methods used in this project draw heavily upon the R package `geomorph` (Adams and Otarola-Castillo 2013).

___

## Introduction

The past few years have witnessed the development of many sophisticated statistical methods that combine geometric morphometrics with phylogenetic comparative methods. In addition, recently developed R packages have made these methods readily accessible to evolutionary biologists. These methods allow researchers to avoid spurious results due to failure to account for the phylogenetic structure of their data (Díaz-Uriarte and Garland 1996; Uyeda et al., 2015), and more excitingly, to test numerous macroevolutionary hypotheses that cannot be addressed with traditional geometric morphometrics or univarite phylogenetic comparative methods (Klingenberg and Marugan-Lobon 2013). 
This document summarizes the analyses and results I've obtained from applying phylogenetic geometric morphometric methods (PGMMs) to a comparative dataset on primate cranial shape from Fleagle et al. (2010). By analyzing this data in a phylogenetic context, I hope to yield insights into the ecological and evolutionary factors underlying variation in primate cranial shape. This investigation is divided into three chapters:

1. Effects of allometry and sexual dimorphism on primate cranial shape evolution
2. Ecomorphology of primate cranial shape
3. Modularity, integration, and rates of primate cranial shape evolution

The same landmark data and phylogeny are used in all three chapters. The primate cranial shape dataset includes 18 landmarks for a male and female specimen from a representative species from most primate genera. Since some species only had one of the sexes sampled, this results in datasets of slightly different sizes (*n* = 64 males, *n* = 61 females). The phylogeny is a consensus tree from 10kTrees (Arnold et al. 2010). The only changes I made to the data involved editing taxa labels so that they match between the data and the phylogeny, and dropping two taxa from the landmark data because they were not present in the phylogeny. Landmarks were aligned using generalized procrustes alignment prior to all analyses, and centroid sizes were recorded for all crania. Here is a view of the species and phylogeny:

![](./figures/SpeciesTree.tiff)

Before doing any phylogenetic analyses of shape, I tested for significant phylogenetic signal using two methods. First, I used a multivariate extension of a well-known univariate measure of phylogenetic signal, Blomberg's K (Adams 2014). Second, I used an approach developed by Klingenberg and Gidaszewski (2010), which uses squared-change parsimony as a test statistic to compare the observed data to a null distribution obtained by permuting data at the tips of the tree. Both methods revealed strong phylogenetic signal in the data (*p* < 0.001), indicating that it is appropriate to analyze this data in a phylogenetic framework. We can visualize the phylogenetic structure of our shape data by projecting the phylogeny and ancestral state reconstructions into tangent space. This is akin to plotting datapoints in the space defined by the first two PC axes, except this plot also depicts the phylogenetic relationships among datapoints and the locations of ancestral states in shape space. 

![](./figures/Phylomorphospace.tiff)

Alright, now that the cranial shape data is ready to go and I have confirmed that the data is phylogenetically structured, its time to ask some more interesting questions. 

## CHAPTER 1: Effects of allometry and sexual dimorphism on primate cranial shape



## CHAPTER 2: Ecomorphology of primate cranial shape



## CHAPTER 3: Modularity, integration, and rates of primate cranial shape evolution

Modularity and integration are closely related concepts that describe patterns of covariation among anatomical structures within organisms. Modular structures are relatively independent of one another, while integrated structures tend to covary. Modularity and integration can be observed at multiple levels: across ontogenetic stages within an individual, across individuals within a population or species, and across species within a higher taxonomic group. In addition, multiple biological processes contribute to patterns of integration at each level of variation, including genetic pleiotropy, linkage disequilibrium, allometry, physical interactions among structures during development, and environmental factors that simultaneously influence the development or evolution of multiple structures. Detecting patterns of modularity and integration and the processes that give rise to them requires different study designs depending on the level of variation and the biological process under investigation.

The vast majority of research on cranial modularity and integration has focused on variation within populations or species (Klingenberg 2013). However, modularity and integration can also be studied on a macroevolutionary scale (Klingenberg & Marugán-Lobón 2013). At this scale, patterns of modularity and integration are influenced by a combination of genetic, developmental, functional, and evolutionary factors. Previous studies have found that overall patterns of cranial modularity are remarkably consistent across primates (bunch of studies, cited in Ackerman 2009) and even across mammals (Goswami 2006; Marroig et al., 2009; Porto et al., 2009). This stability is likely due in part to genetic and developmental constraints, and in part due to stabilizing selection. On the other hand, while *patterns* of modularity are largely stable across taxa, the *magnitude* of integration across modules can vary considerably (e.g., Claverie & Patek, 2013), and it has been hypothesized that variation in the strength of integration can influence morphological evolution on a macroevolutionary scale (Olson & Miller, 1958). For instance, it has been hypothesized that modularity will tend to increase through evolutionary time as natural selection favors the relaxation of constraints that integration places on morphological diversification (Marcot & McShea, 2007; Porto et al., 2009). Another hypothesis is that greater modularity facilitates higher rates of evolutionary change, because different modules can respond independently to different selection pressures, while greater integration results in slower evolutionary rates (Goswami & Polly, 2010; Claverie & Patek, 2013; Goswami et al., 2015). 




## Conclusions

This project demonstrates the broad scope of macroevolutionary questions that can be addressed when geometric morphometrics and phylogenetic comparative methods are integrated. 

## References

- Adams DC. 2014. A generalized K statistic for estimating phylogenetic signal from shape and other high-dimensional multivariate data. Syst Biol 63:685-697.

- Adams DC, Otarola-Castillo E. 2013. geomorph: an R package for the collection and analysis of geometric morphometric shape data. Meth Ecol Evol 4:393-399.

- Arnold C, Matthews LJ, Nunn CL. 2010. The 10kTrees Website: a new online resource for primate phylogeny, v3. Evol Anthropol 19: 114–118.

- Díaz-Uriarte R, Garland T. 1996. Testing hypotheses of correlated evolution using phylogenetically independent contrasts: sensitivity to deviations from Brownian Motion. Syst Biol 45(1):27-47.

- Fleagle JG, Gilbert CC, Baden AL. 2010. Primate cranial diversity. Am J Phys Anth 142:565-578.

- Klingenberg CP, Gidaszewski NA. 2010. Testing and quantifying phylogenetic signals and homoplasy in morphometric data. Syst Biol 59:245-261.

- Klingenberg CP, Marugan-Lobon J. 2013. Evolutionary covariation in geometric morphometric data: analyzing integration, modularity, and allometry in a phylogenetic context. Syst Biol 62(4):591-610.

- Uyeda JC, Caetano DS, Pennell MW. 2015. Comparative analyses of principal components can be misleading. Syst Biol 64(4):677-89.

___