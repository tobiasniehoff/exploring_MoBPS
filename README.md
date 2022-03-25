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
The excercises in this section are partially my own or from the MoBPS workshop by Torsten Pook. Code, solutions and comments are my own.

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
