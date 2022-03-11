# exploring_MoBPS
Example scipts and additional functions using [MoBPS](https://github.com/tpook92/MoBPS).
This respository is meant to share scripts and functions that have been built with and/or for MoBPS. Contributions will mostly be related to animal breeding. I will upload something every now and then from my own projects. I do not guarantee that code will work as described or intended and always recommend testing it for your own purposes.

Please et me know if you spot mistakes, have questions or just want to say Hi! :)

## Scripts
This is a collection of scripts of things you can do with MoBPS.

`simulate_population_history_Jibrila_et_al_2020.R`
This script simulates a historic population as described in [Jibrila et a. \(2020\)](https://doi.org/10.1186/s12711-020-00562-6). Comments are made in the script. In brief, a random mating population is simulated over 3000 generations with changes in population size. This script could run on my laptop and needed some hours to finish.

## Functions
This folder contains functions made for and/or with MoBPS. Some functionality may be incorporated in future MoBPS versions. The fuctions are described more in the file. Examples are given at the end of function definitions.

**`min.coanc.mating`**
- This function does minimum coancestry mating. It is a wrapper around the function `matings()` of the package [`optiSel`](https://doi.org/10.1186/s12859-018-2450-5). It uses pedigree relationship and assumes equal contributions within sex. Using `Rsymphony::Rsymphony_solve_LP` is way faster than using the default solver and is less prone to fail. Installing `RSymphony` on the cluster is a nightmare. I found that using `lpsymphony` is easier to [install](https://www.bioconductor.org/packages/release/bioc/html/lpsymphony.html) and does the same.
If no solver works or you don't want to wait, my personal suggestion is to use the `avoid.mating.fullsib` and `avoid.mating.halfsib` functionality in MoBPS as most inbreeding and reduction in within family variance can be expected from such crosses. I think this is discussed in [Sonesson & Meuwissen (2000)](https://doi.org/10.1186/1297-9686-32-3-231).

**`get.direct.ancestors`**
- Gives you a pedigree containing all ancestors up to a given numbe rof generations back. 
