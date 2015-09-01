# Primate cranial shape evolution: Combining geometric morphometrics and phylogenetic comparative methods

[Randi H. Griffin]()

This repository contains the R code for my ongoing dissertation work. 

___

### Data

The landmark data come from a study published in 2010 by J.G. Fleagle, C.C. Gilbert, and A.L. Baden (Fleagle et al., 2010). Hopefully I will get permission to add that file soon.

`PrimateTree.nex` is a consensus phylogeny downloaded from Version 3 of 10kTrees (Arnold et al., 2010), pruned to taxa present in the landmark data.

`BodyMass.csv` contains male and female body mass estimates from Smith and Jungers (1994). 

`EcologyData2.csv` and`EcologyData3.csv` contain character states for activity pattern, diet, and locomotion style for the primates in our analysis. The two datasets code activity pattern and locomotion styles as 2 and 3 state characters, respectively. These data were obtained from the following public databases:

### Analysis

`Analysis.R` contains the code to run the analyses. Unfortunately, the code is useless without the cranial landmark data.

### Results

`Synopsis.md` contains an overview of my analyses and results so far. 

The folder **figures** contains the figures that are included in `Synopsis.md`. 

### References

- Arnold C., Matthews L.J., and C.L. Nunn. The 10kTrees Website: a new online resource for primate phylogeny, v3. Evol Anthropol. 2010;19: 114–118.

- Fleagle, J.G., Gilbert, C.C., and A.L. Baden. 2010. Primate cranial diversity. Am. J. Phys.
Anth. 142:565-578.

___