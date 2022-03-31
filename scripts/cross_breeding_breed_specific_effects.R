# MoBPS_1.8.09

library(MoBPS)

# reale haplotypen
# row corresponds to a SNP
# column corresponds to a haplotype of an individual
# So each diploid individual will have two columns
data1 <- matrix(rbinom(80,1, 0.1245), nrow=10)
data2 <- matrix(rbinom(60,1, 0.7), nrow=10)

# reale map
map <- cbind(1, paste0("SNP", 1:10), seq(0.1,1, length.out=10), 1:10*100)
colnames(map) <- c("chr", "SNP name", "Morgan position", "bp")

map <- rbind(map, map, map)

map[1:10,2] <- paste0(map[1:10,2], "_BreedA")
map[1:10+10,2] <- paste0(map[1:10+10,2], "_BreedB")
map[1:10+20,2] <- paste0(map[1:10+20,2], "_Array")

# zusammengefasster datensatz
dataset <- matrix(0, nrow=30, ncol=14)
dataset[1:10,1:8] <- data1
dataset[11:20,9:14] <- data2
dataset[21:30,1:14] <- cbind(data1,data2)


# Adding breedA to the population object
population <- creating.diploid(dataset = dataset[,1:8], map = map,
                               sex.s = c(1,1,2,2), name.cohort = "BreedA_1")
get.geno(population, cohorts = c("BreedA_1_M", "BreedA_1_F"))
get.map(population)

# Adding breedB to the population object. Animals of BreedB should be in the same generation
population <- creating.diploid(population = population,
                               dataset = dataset[,9:14], map = map,
                               sex.s = c(1,2,2), name.cohort = "BreedB_1"
)
get.geno(population, cohorts = c("BreedA_1_M", "BreedA_1_F"))
get.geno(population, cohorts = c("BreedB_1_M", "BreedB_1_F"))

get.geno(population, gen = 1)

# showing markers in order of added
reordered <- unlist(lapply(1:3, FUN = function(number)seq(from=number,to=30,by=3)))
get.geno(population, gen = 1)[reordered,]

population <- add.array(population, 
                        marker.included = rep(c(FALSE, FALSE, TRUE),10),
                        array.name = "All markers")

population$info$snp.name
population$info$snp.position
population$info$array.name

# now adding marker effects
# the whole point of this allele-fixing approach is that the user can upload 
# marker effects that are breed specific
#
# let's say we estimated effects for markers in both real breeds separately

# These are the effects we use for BreedA
# for now, let's don't care about the effect sizes, additivity, dominance or epistasis
# this is just an example
# I will here assume perfect additivity
#
# Also, maybe we don't have effect estimates for all SNPs
# This can happen if an allele is fixed in the population or if the locus
# wasn't included in the effect size estimation but should be simulated
# I will here set the effect size of the last SNP in both Breeds to 0 to mimick this

SNP <- 1:10
chromosome <- 1
effect0 <- rep(0, times = 10)
effect1 <- c(abs(rnorm(9, mean=0, sd=5)),0)
effect2 <- effect1*2
effect.matrix.BreedA <- cbind(SNP, chromosome, effect0, effect1, effect2)

# for testing, we will set the effect sizes of markers of animals
# from BreedB to much larger values
SNP <- 1:10
chromosome <- 1
effect0 <- rep(0, times = 10)
effect1 <- c(abs(rnorm(9, mean=100, sd=5)),0)
effect2 <- effect1*2
effect.matrix.BreedB <- cbind(SNP, chromosome, effect0, effect1, effect2)

effect.matrix.Both <- cbind(SNP, chromosome, matrix(0, nrow = 10, ncol = 3))

effect.matrix <- rbind(effect.matrix.BreedA, effect.matrix.BreedB, effect.matrix.Both)

# as can be seen here, we need to simulate 30 markers
nrow(get.geno(population, gen = 1))

# reordering the effect matrix to match the positions
# this is the order the markers should have
get.geno(population, gen = 1)
v <- c()
for(i in 0:9){v <- c(v, c(1,11,21)+i)}

effect.matrix <- effect.matrix[v,]

# we need to change the SNP position in the effect matrix
effect.matrix[,1] <- 1:30

pop_store <- population
population <- pop_store

# important to set base.bv=0. Otherwise the default of 100 is used which does 
# not make sense if effects are specified. Effectively, this would add 10 to 
# the sum of the allele effects of every individual
population <- creating.trait(population, real.bv.add = effect.matrix, base.bv = 0)
get.bv(population, gen = 1)
get.qtl(population)
get.qtl.effects(population)

###
# Here we are checking if the correct tbv is reported
geno1 <- get.geno(population, cohorts = c("BreedA_1_M", "BreedA_1_F"), use.id = T)[reordered,][1:10,]
# this is our own calculation of marker effects multiplied with allele count
sum(geno1[,1]*effect.matrix.BreedA[,4])
# this is the MoBPS output on true breeding values
# both should be identical
get.bv(population, database = get.database(population, id = colnames(geno1)[1]))

# checking with QTL effects stored in the population object
sum(geno1[,1]*get.qtl.effects(population)[[1]][[1]][seq(from=1,to=30,by=3),4])
###

################################################################################
# Now trying some breeding and mixing of the populations

population <- breeding.diploid(population, breeding.size = 100, 
                               selection.f.cohorts = "BreedA_1_F",
                               selection.m.cohorts = "BreedA_1_M", 
                               name.cohort = "BreedA_2")
get.cohorts(population)
get.geno(population, gen = 2)

for(i in 2:5){
  population <- breeding.diploid(population, breeding.size = 100, 
                                 selection.f.cohorts = paste0("BreedA_", i, "_F"),
                                 selection.m.cohorts = paste0("BreedA_", i, "_M"), 
                                 selection.criteria = c("bv", "bv"), 
                                 name.cohort = paste0("BreedA_", i+1), 
                                 selection.size = c(10, 20))
}

get.bv(population, cohorts = c("BreedA_1_M", "BreedA_1_F"))
get.bv(population, gen = length(population$breeding))

# AF in the first animals of BreedA
rowMeans(get.geno(population, cohorts = c("BreedA_1_M", "BreedA_1_F"))[reordered,][1:10,])/2
# AF in the last generation of BreedA
rowMeans(get.geno(population, gen = length(population$breeding))[reordered,][1:10,])/2

get.cohorts(population)

# here, hybrids between BreedA and BreedB are created
population <- breeding.diploid(population, 
                 selection.f.cohorts = "BreedB_1_F",
                 selection.m.cohorts = "BreedA_6_M", breeding.size = 200, 
                 name.cohort = "F1")

get.bv(population, gen = length(population$breeding))

# value for sigma.e is a random choice
# Using a heritability would be confusing as the breed differences inflate
# the genetic variance and thus also the environmental variance if a 
# target heritability is to be met
population <- breeding.diploid(population, 
                               phenotyping.gen = length(population$breeding), 
                               sigma.e = 5)
# Phenotypes of hybrids can be derived
get.pheno(population, gen = length(population$breeding))


# creating some interbred animals
for(i in 1:4){
  population <- breeding.diploid(population, breeding.size = 200)
}

population <- breeding.diploid(population, 
                               phenotyping.gen = length(population$breeding), 
                               sigma.e = 50)
# Phenotypes of crossbreds can be derived
get.pheno(population, gen = length(population$breeding))
# you can see that some animals have negative phenotypes
# this is because BreedA has a mean of ~15 and BreedB has a mean of 1000.
# If intercrossed, some animals are genetically closer to BreedA by chance 
# and thus also have a lower phenotype - by the environmental standard deviation
# is huge in comparison


# this is how you would do breeding value estimation
# you would have to use the "All markers" array
# here, some purebred aniamls from BreedA, hybrids and crossbreds are used
population <- breeding.diploid(population, bve = TRUE, 
                               bve.array = "All markers", 
                               bve.gen = 5:9, 
                               relationship.matrix = "vanRaden")

get.bve(population, gen = 5:9)

### Bonus
# with the function get.breed.proportion() we can retrieve the breed proportion 
# of every animal
# animals in group 2 must be founders
# the function is defined in MoBPS_exploration_functions.R
source('../functions/MoBPS_exploration_functions.R')

breed.prop <- get.breed.proportion(population,
                     gen1 = length(population$breeding),
                     cohorts2 = c("BreedA_1_M", "BreedA_1_F"))

# the mean proportion of the genome that can be traced back to 
# breedA animals is about 50%
mean(breed.prop)


### 
# Here, we will simulate some generations of phenotypic selection
population <- breeding.diploid(population, 
                               phenotyping.gen = length(population$breeding), 
                               sigma.e = 50)

# this is the phenotypic mean before selection
(pheno.mean.before <- mean(get.pheno(population, gen = length(population$breeding))))
# this is the mean of the proportion of breedA genome in the population before selection
(breed.prop.before <- mean(breed.prop))

for(i in 1:5){
  # here a new generation is generated
  population <- breeding.diploid(population, breeding.size = 200, 
                                 selection.size = c(10,25),
                                 selection.criteria = c("pheno", "pheno"), 
                                 name.cohort = "Cross")
  
  # here pehnotyping is done
  population <- breeding.diploid(population, 
                                 phenotyping.gen = length(population$breeding), 
                                 sigma.e = 50)
}

breed.prop.cross <- get.breed.proportion(population,
                     gen1 = length(population$breeding),
                     cohorts2 = c("BreedA_1_M", "BreedA_1_F"))


# this is the phenotypic mean before selection
(pheno.mean.after <- mean(get.pheno(population, gen = length(population$breeding))))
# this is the mean of the proportion of breedA genome in the population before selection
# you can see that the contribution of BreedA to the current generation has dropped dramatically
# this is because the allele effects of BreedA are inferior compared to BreedB
(breed.prop.after <- mean(breed.prop.cross))
