# Why this repository?
I am working with MoBPS all the time and wanted to have a place to share code and ideas. It can take a while to set up your idea. I hope that this repository can help or give inspiration especially when you are just starting working with MoBPS.

<!-- 
# My motivation:
I am interested in having as many people as possible using MoBPS. This ensures continued development which benefits all users. With this repository, I want to help help increasing the user base. I see MoBPS as the best tool for breeding program simulation for my projects. 
Please reach out if you want to have a chat.
-->

# What is this repository?
Example scipts and additional functions using [MoBPS](https://github.com/tpook92/MoBPS).
This respository is meant to share scripts and functions that have been built with and/or for MoBPS. Contributions will mostly be related to animal breeding. I will upload something every now and then from my own projects. I do not guarantee that code will work as described or intended and always recommend testing it for your own purposes.

There are 3 sections in this repo:
- **Exercises**: Here you can find exercises for practice. Th ecorresponding folder contains my solutions.
- **Scripts**: Here you can find scripts of my previous projects or examples for some implementations.
- **Functions**: Here, only the folder is relevant. You can find functions that extend the functionality of MoBPS.
<!--
- **General remarks**: This section does not contain code. I am commenting on open questions I have or give general recommendations.
-->
<!-- 
Please contact me know if you spot mistakes, have questions or just want to say Hi! :)
-->

## Version
- MoBPS actions can be sped up with the packages `miraculix` [https://rdrr.io/cran/miraculix/src/R/auxiliary.R](https://rdrr.io/cran/miraculix/src/R/auxiliary.R) and `RandomFieldUtils` [https://cran.r-project.org/web/packages/RandomFieldsUtils/index.html](https://cran.r-project.org/web/packages/RandomFieldsUtils/index.html). Certain versions of the packages do not work on Windoes or Linux. This might change in the future as the packages are actively developed. For now, the most recent stable versions for <br/>
 **Windows**:<br/>
 `RandomFieldsUtils_0.6.6`<br/>
 `miraculix_1.0.0.1`<br/>
 **Linux**:<br/>
 `RandomFieldsUtils_1.0.6`<br/>
 `miraculix_1.0.5`
 
 These versions can be downloaded on Torsten's MoBPS GitHub page [https://github.com/tpook92/MoBPS](https://github.com/tpook92/MoBPS).
 
 ```
Installation on Linux:
install.packages("MoBPS_1.8.07.tar.gz", repos = NULL, type = "source")
install.packages("RandomFieldsUtils_1.0.6.tar.gz", repos = NULL, type = "source")
install.packages("miraculix_1.0.5.tar.gz", repos = NULL, type = "source")
```
**NOTE** MoBPS is continously updated. Thus, some scripts may not work with a certain version of MoBPS. I try my best to always mention the version with which a script was tested in the `.R` files.

**NOTE** MoBPS, RandomFieldsUtils and miraculic do **NOT** work in R versions **`R4.3.1+`** (status 26.06.2024). To my knowledge, Torsten Pook and Martin Schlather (the developer of RandomFieldUtils and miraculix) is working on this.

## Exercises
The excercises in this section are partially my own or from the MoBPS workshop by Torsten Pook. Code, solutions and comments are my own. Note that there is not one solution only to implement a program in MoBPS.

### Hybrid breeding
You will learn how to read in real marker data from vcf files, create traits and add restrictions to matings. The file with solutions is called **`task6_solution_hybrid_breeding.R`**. This excerise is adapted from the MoBPS workshop (March 2022). 
- Generate a founder population with 10 individuals from two different gene pools
- Genotypic data for the individuals is given via Pool1.vcf and Pool2.vcf
- Add three traits:<br/>
- One with 10 purely additive QTLs<br/>
– One with 1000 purely additive QTLs<br/>
– One with 500 purely additive QTLs and 500 dominant QTLs of equal size<br/>
– Traits should be uncorrelated<br/>
– For all traits phenotypic mean for the founders should be 100 with a variance of 5<br/>
- Generate 100 offspring by random mating between individuals from the two gene pools
- Compare the genomic values of the parents and offspring
- How often was each individual used for reproduction?
- Generate a PCA for all simulated individuals

### Breeding value estimation
You will learn how to do PBLUP, GBLUP and ssGBLUP in MoBPS and how to use different marker arrays in bve. The file with solutions is called **`task7_solution_BreedingValueEstimation.R`**. This excerise is adapted from the MoBPS workshop (March 2022). This excercise  is split into two subparts. If you are struggling with the first part you can load in the `.Rdata` object with an already generated population<br/>
- Generate a population list with 12 generations:<br/>
    - Each generation contains 50 males, 50 female with parents of the previous generation<br/>
    - Use a genetic map with 25.000 SNPs, 5 chromosomes with a length of 3 Morgan each
    - Generate a trait with heritability of 0.3 and 1‘000 purely additive QTLs
    - Make sure that:
        - In generation 10 only males are phenotyped
        - In generation 11 & 12 all individuals are phenotyped
        - In generation 10 & 11 all individuals are genotyped
        - In generation 12 only 20% of all individuals are genotyped
- Perform breeding value estimations for individuals in generation 10, 11, 12 or combinations of the cohorts. Use:
    - Genomic breeding value estimation
    - Pedigree based breeding value estimation
    - Single step
    - Assume individuals are only genotyped for 10‘000 / 100 randomly selected markers
- Generate a plot to showcase real genomic values and breeding values for generation 10

### Index selection
You will use how to do selection with an index for multiple traits. Index selection is done in most contemporary breeding programs. In this exercise, selection should be done for an arbitrary index comprised of 3 correlated traits. The file with solutions is called **`index_selection.R`**.
- Simulate a population with 5000 SNPs on 28 chromosomes and 300 individuals
- Create 3 traits and make up some genetic correlations. Create the correlation matrix.
- Make up some genetic variances and heritabilities. Based on both, calculate error variance.
- Add the three traits to the population. Each trait should have 1000 causla additive loci.
- Compare your input correlation matrix with the correlation in the population.
- Let the population mate randomly for some generations to build up some LD.
- Do a PBLUP breeding value estimation based on phenotypes from only females.
- Do index selection for all animals in the last generation based on ebvs. Pay attention to the variance you standardize ebvs with.

### Offspring phenotypes
You will learn how to use offspring phenotypes. This is useful in cases the animal of interest (e.g. a rooster) does not have a phenotype (e.g. total egg mass). The file with solutions is called **`task8_solution_OffpringPhenotypes.R`**. This excerise is adapted from the MoBPS workshop (March 2022).
- Simulate a population with 10 male individuals and 90 female individuals
- Generate a single trait with 1000 purely additive QTL
- Generate 45 male and 45 female offspring
- Generate phenotypes for all female offspring
- Calculate the average phenotype of the offspring from each respective founder
- As comparison perform a breeding value estimation and compare which of the two selection criteria is more suited to estimate the underlying true genomic value

## Scripts
This is a collection of scripts of things you can do with MoBPS.

**`short_introduction_MoBPS.R`**
- This script gives a short introduction into MoBPS and shows aspects one is often interested in such as correlated traits, different breeding value estimation methods and a simple breeding scheme with overlapping generations. It is meant for people just starting with MoBPS.<br/>
![Age_proportion_parents](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/Age_proportion_parents.png)

**`cross_breeding_breed_specifc_effects.R`**
- In this script, I explore how one can simulate crossbreeding when the known effect sizes for alleles differ per breed/population. This is the case when they are estimated separately per breed (i.e., the effect size of allele A in breed A is different from the effect size of allele A in breed B) and one wants to use the "real" estimated allele effects as input for simulation. The way this is implemented here, is by placing a SNP for breed A on the same position as for breed B. Note that for this only the genetic position (cM), not the physical position (bp) matters. Thus, an animal can only have one allele at a locus - either the one from breed A or from breed B because no recombination can happend between them. Then, the allele frequency for the SNP from breed B is set to 0 in breed A animals and vice versa. In breeding value evalutations, an array conataining both SNPs is needed. For other genomic evaluations, make sure you use the appropriate marker array.<br/>
Note that MoBPS does not have a built in strategy to consider the breed-origin-of-alleles in breeding value estimation. For an extensive review on considering crossbreed data in genomic prediction (also covering the inclusion of the breed-of-origin for alleles), please see [Duenk et al. \(2021\)](https://doi.org/10.1093/jas/skab205).<br/>
The scripts further explores how selection changes breed proportions in a crossbred population.<br/>
This method is complicated. Compatibility with other MoBPS functions is not guaranteed.

**`connectedness_population.R`**
- This script simulates a breeding program that is having two separate herds. Some semen is exchanged between the herds every generation. This is implemented as such: e.g. 10% of all father for herd A are from herd B. These sires are the very best ones in herd B. And the same is done for sires of herd B. The script calculates and plots the average true breeding value and the Fst value. I have not seen a huge population differentiation during development.

**`simulate_population_history_Jibrila_et_al_2020.R`**
- This script simulates a historic population as described in [Jibrila et al. \(2020\)](https://doi.org/10.1186/s12711-020-00562-6). Comments are made in the script. In brief, a random mating population is simulated over 3000 generations with changes in population size. This script could run on my laptop and needed some hours to finish.

**Optimum Genetic Contribution considering the number of to be selected parents**
- The code and scripts for this investigation are in the folder `OGC_number_of_parents_script`. I compare different strategies using optimum contributions as a criterion to make selection decisions when only a limited number of parents can be selected for the next generation. I am using [optiSel](https://cran.r-project.org/web/packages/optiSel/index.html). The number of parents is considered by calculating OC for all selection candidates and then remove the animal with the lowest OC. Then, optimum contributions are recalculated and again the animal with the lowest OC is removed. This is repeated until only the desired number of animals remains with non-zero contributions. I implemented this is a wrapper function. <br/>
Since this one-by-one removal can take some time, I also wrote an algorithm that significantly speeds up the process by removing more than one animal at a time. The script uses the excel sheet as an input file. You can change the nuber of selected parents, whether EBV should be maximised or kinship minimized, a minimum genetic gain or maximum inbreeding increase (by specifying ne) and more. OGC styles are 1) selection based on BV and contributions of selected animals are optimized, 2) optimum contributions are calculated for all animals and animals with highest contributions are selected and 3) the iterative approach of recalculating contributions and removing one animal every iteration.<br/>
An `.sh` is provided so that analysis could be run on a cluster.<br/>
Some findings based on this script were presented on 14.10.2022 at the CiBreed conference in Göttingen.<br/>
`Note:` tools like [AlphaMate](https://alphagenes.roslin.ed.ac.uk/wp/software-2/alphamate/) or [mateSel](https://www.matesel.com/) use an evolutionary algorithm to find the best solution and can consider the number of to be selected parents as a constraint directly.

![inbreeding_OGC_strategies](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/ptab9%201251-1500%20wo%202.2%20inbreeding.png)

**Integrating physical marker postions (bp) into genetic positions (Morgan)**
- The script is in the folder `Integrating genetic map`. The genetic map needs to be downloaded from supplemenatry material of [Groenen et al. \(2009\)](10.1101/gr.086538.108): A high-density SNP-based linkage map of the chicken genome reveals sequence features correlated with recombination rate.
- I used MoBPS version 1.11.40
- This script contains functionality to derive genetic positions (in Morgan) for markers for which only physical positions (in bp) are available. In this example, I am using the genetic map estimated for chicken by [Groenen et al. \(2009\)](10.1101/gr.086538.108). That map contains 8621 markers on autosomes. I want to derive genetic positions for 600K SNPs that are on the Affymetrix Chicken600K Array (can be found in the packge MoBPSmaps `MoBPSmaps::map_chicken1`). The MoBPS manual recommends a conversion factor of 30.000.000 bp per Morgan for chicken but this not reflecting chromosome or region specific differences in recombination rate. That conversion can differ between chromosomes can be visualized by plotting the physical length (bp) per genetic length (M) of every chromosome based on teh daa provided by [Groenen et al. \(2009\)](10.1101/gr.086538.108) (chromosome 16 and 25 have no height because not enough SNPs were left after quality checking to estimate length reliably).<br/>
![chromosome specific vonversion factor](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/scripts/Integrating_genetic_map/chromosome_specific_conversion_fatcor.png)

- So, I would like to be as accurate as possible by estimating a conversion function based on the genetic map by [Groenen et al. \(2009\)](10.1101/gr.086538.108). This works by estimating a linear line between every pair of markers with function `stats::splinefun()`, method `hyman`. The plot below shows such a conversion function for chromosome 15.<br/>
![estimated spline function](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/scripts/Integrating_genetic_map/positions_groenen.png)

- Once this conversion function is estimated, the physical positions can be provided as input and the genetic position will be assigned based on the estimated function.<br/>
![added_600K_SNPs](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/scripts/Integrating_genetic_map/added_600K_SNPs.png)

- The function `function.bp.Morgan()` (provided in the script) also does some quality control. Duplicated SNPs are removed, SNPs with rank issues are removed (bp rank higher than M rank) and outliers are removed. For example, there are a couple of SNPs on chromosome 22 whose Morgan position does not match with the Morgan position of SNPs that are neighbors on the bp scale (plot below). This could be because the ranking (genome build) of the physical map is incorrect. Here, I remove them by estimating a linear regression line and removing points that deviate by more than 4 OR the SD that corresponds to 1/(number SNPS chromosome quantile standard deviations (whichever is larger) from the prediction. For chromosome 22, this means that the outlier SNPs at around 0.4 M are removed.<br/>
![added_600K_SNPs](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/scripts/Integrating_genetic_map/chrom22_raw.png)

- Functionality to integrate SNPs that are beyond the range of the template genetic map is also implemented as well as functionality that assigns genetic positions purely based on a predefined conversion factor.

## Functions
This folder contains functions made for and/or with MoBPS. Some functionality may be incorporated in future MoBPS versions. The fuctions are described more in the file. Examples are given at the end of function definitions.<br/>
**Functionality for the Optimum Genetic Contribution script can be found under Scripts -> OGC_number_of_parents_script -> `functions_for_OGC_script.R`**


**`get.direct.ancestors()`**
- Gives you a pedigree containing all ancestors up to a given number of generations back. 

**`min.coanc.mating()`**
- This function does minimum coancestry mating. It is a wrapper around the function `matings()` of the package [`optiSel`](https://doi.org/10.1186/s12859-018-2450-5). It uses pedigree relationship and assumes equal contributions within sex. Using `Rsymphony::Rsymphony_solve_LP` is way faster than using the default solver and is less prone to fail. Installing `RSymphony` on the cluster is a nightmare. I found that using `lpsymphony` is easier to [install](https://www.bioconductor.org/packages/release/bioc/html/lpsymphony.html) and does the same.
If no solver works or you don't want to wait, my personal suggestion is to use the `avoid.mating.fullsib` and `avoid.mating.halfsib` functionality in MoBPS as most inbreeding and reduction in within family variance can be expected from such crosses. I think this is discussed in [Sonesson & Meuwissen (2000)](https://doi.org/10.1186/1297-9686-32-3-231).

**`get.full.database()`**
- This function returns a database with every animal on a single row, i.e., the number of rows is identical to the number of animal in the database. This database makes some handeling easier. It is acceptable by all MoBPS functions. I have not tested the case where clones or copies are used but I suspect that it would work the same.

**`get.parental.average()`**
- This function returns the parental average of every individual. The average of phenotypesm estimated breeding values and true breeding values can be calculated by specifying "pheno", "ebv" or "bv". I have only tested this function for the case that one trait is used and only one phenotype record of that tait is available per animal. **Attention plant breeders:** the output will not be meaningful for you if the plant of interest are derived by selfing. If you have selfed plants, you probably want to make sure that the original lines that were used to create the F1 are use to derive th parental average. This is ***not*** yet implemented.

**`get.corrected.bv()`**
- This function returns the corrected bvs of animals. Options to correct by are "generation.mean" and "parental.mean". For some cases, the genetic variance based on corrected breeding values (here: true breeding values) is of interest. I have not tested this function for multiple traits. internally, it calls `get.parental.average()`.

**`get.corrected.bve()`**
- This function works the same as `get.corrected.bv()` - just for estimated breeding values.

**`get.breed.proportion()`**
- This function checks what proportion of the genome of an animal is derived from a specific group of founders. If the group of founder is chosen to contain all but only animals from a certain breed/poplation, the value given as output is the breed proportion. The function only works if the animals in group 2 contain only founder animals.

<!-- 
https://gist.github.com/jonikarppinen/47dc8c1d7ab7e911f4c9
## General remarks
-Simulation tools, also MoBPS, offer the user to input heritabilities which is then used together with either the environmental or the genetic variance to derive the genetic or the environmental variance, respectively. This is coming from a way of thinking that that a heritability is a property of a trait. However, a heritability is a property of the environment AND the population. This is important in simulations when following a population over generations. Assume your breeding program is constantly increasing a trait. Assume the heritability was estimated based on a single generation and this is the heritability you want to use in simulation or e.g. phenotypes. 
There are two problems:
1) When simulating phenotypes, a tool (including MoBPS) would derive the true genetic variance (which is known to the tool). Together with the heritability, the environmental variance will be estimated. The environmental variance is then used to draw random error effects. When simulating many generations, the genetic variance can be expected to decrease. In a real population, we would also assume the heritability to decrease. With the approach explained here, the simulation tool would however also decrease the environmental variance. This can be avoided in MoBPS by providing the environmental standard deviation instead of the heritbaility when phenotyping. It should be noted though, that unless you use ridiculously small population sizes, the genetic variance probably does not change much in only a few generations.
2) Again, assume the used heritability is estimated based on a single generation. Now, if you want to phenotype several generations at once, the simulation tool will overestimate the genetic variance and in turn also the environmental variance. This is because the generations have different means for the trait. To avoid this problem, one can use the parameter `variance.correction` in MoBPS. In most cases, it is probably bets to chose `generation.mean` instead of `parental.mean`. The latter corrects "too much" and removes the part of the variance that is due to the mating strategy. Alternatively, one could provide the environmental variance directly for drawing random effects. The problem is visualized in the figure ![variance_problem](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/problem_variances.png).

Sidenote: I am currently pondering about environmental variances. Imagine this: Your population has quite a low mean, say 10. The environmental standard deviation in this population is, say 4. After many generations, the trait was increased to 10,000. Would the population then still have an environmental sd of 4 or would it increase? 
The question behind this is whether the variance is relative to the tbv or not. If you know the answer or have experience, please contact me.
-->
