# Why this repository?
I am working with MoBPS all the time and wanted to have a place to share code and ideas. It can take a while to set up your idea. I hope that this repository can help or give inspiration especially when you are just starting working with MoBPS.

# My motivation:
I am interested in having as many people as possible using MoBPS. This ensures continued development which benefits all users. With this repository, I want to help help increasing the user base. I see MoBPS as the best tool for breeding program simulation for my projects. 
Please reach out if you want to have a chat.

# What is this repository?
Example scipts and additional functions using [MoBPS](https://github.com/tpook92/MoBPS).
This respository is meant to share scripts and functions that have been built with and/or for MoBPS. Contributions will mostly be related to animal breeding. I will upload something every now and then from my own projects. I do not guarantee that code will work as described or intended and always recommend testing it for your own purposes.

Please contact me know if you spot mistakes, have questions or just want to say Hi! :)

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

**`connectedness_population.R`**
- This script simulates a breeding program that is having two separate herds. Some semen is exchanged between the herds every generation. This is implemented as such: e.g. 10% of all father for herd A are from herd B. These sires are the very best ones in herd B. And the same is done for sires of herd B. The script calculates and plots the average true breeding value and the Fst value. I have not seen a huge population differentiation during development.

**`simulate_population_history_Jibrila_et_al_2020.R`**
- This script simulates a historic population as described in [Jibrila et al. \(2020\)](https://doi.org/10.1186/s12711-020-00562-6). Comments are made in the script. In brief, a random mating population is simulated over 3000 generations with changes in population size. This script could run on my laptop and needed some hours to finish.


## Functions
This folder contains functions made for and/or with MoBPS. Some functionality may be incorporated in future MoBPS versions. The fuctions are described more in the file. Examples are given at the end of function definitions.

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

## General remarks
-Simulation tools, also MoBPS, offer the user to input heritabilities which is then used together with either the environmental or the genetic variance to derive the genetic or the environmental variance, respectively. This is coming from a way of thinking that that a heritability is a property of a trait. However, a heritability is a property of the environment AND the population. This is important in simulations when following a population over generations. Assume your breeding program is constantly increasing a trait. Assume the heritability was estimated based on a single generation and this is the heritability you want to use in simulation or e.g. phenotypes. 
There are two problems:
1) When simulating phenotypes, a tool (including MoBPS) would derive the true genetic variance (which is known to the tool). Together with the heritability, the environmental variance will be estimated. The environmental variance is then used to draw random error effects. When simulating many generations, the genetic variance can be expected to decrease. In a real population, we would also assume the heritability to decrease. With the approach explained here, the simulation tool would however also decrease the environmental variance. This can be avoided in MoBPS by providing the environmental standard deviation instead of the heritbaility when phenotyping. It should be noted though, that unless you use ridiculously small population sizes, the genetic variance probably does not change much in only a few generations.
2) Again, assume the used heritability is estimated based on a single generation. Now, if you want to phenotype several generations at once, the simulation tool will overestimate the genetic variance and in turn also the environmental variance. This is because the generations have different means for the trait. To avoid this problem, one can use the parameter `variance.correction` in MoBPS. In most cases, it is probably bets to chose `generation.mean` instead of `parental.mean`. The latter corrects "too much" and removes the part of the variance that is due to the mating strategy. Alternatively, one could provide the environmental variance directly for drawing random effects. The problem is visualized in the figure ![variance_problem](https://github.com/tobiasniehoff/exploring_MoBPS/blob/main/problem_variances.png).

Sidenote: I am currently pondering about environmental variances. Imagine this: Your population has quite a low mean, say 10. The environmental standard deviation in this population is, say 4. After many generations, the trait was increased to 10,000. Would the population then still have an environmental sd of 4 or would it increase? 
The question behind this is whether the variance is relative to the tbv or not. If you know the answer or have experience, please contact me.
