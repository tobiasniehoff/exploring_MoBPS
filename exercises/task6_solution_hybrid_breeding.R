# R version 4.0.2 (2020-06-22)
# MoBPS_1.8.07

library(MoBPS)

pop <- creating.diploid(vcf = "Pool1.vcf", name.cohort = "Pool1_Founder",
                        bpcm.conversion = 1000000)
get.database(pop, gen = 1)

pop <- creating.diploid(population = pop, vcf = "Pool2.vcf", 
                        name.cohort = "Pool2_Founder",
                        bpcm.conversion = 1000000)
# bpcm.conversion = 1000000 is for converting base pairs to (centi)Morgan
# as vcf files do only store the bp position
# the conversion could be checked with get.map()
# for MoBPS version 1.8.07 however, get.map() does not give meaningful positions
# for genetic positions
# 1 Morgan = 100.000.000 bp is a common conversion rate  (e.g. is assumed in BEAGLE)
# This can differ between species (e.g. chicken: ~30.000.000 bp) 
# or be dependent between telomere / centromere

summary(pop)
get.database(pop, gen = 1)
get.cohorts(pop)
pop$info$chromosome # the species has 5 chromosomes

length(pop$breeding) # numbers of generations in the pop object

# Adding traits with 10 purely additive QTLs
pop <- creating.trait(pop, n.additive = 10, 
                      trait.name = "AdditiveSmallTrait",
                      mean.target = 100, var.target = 5)

# Adding traits with 1000 purely additive QTLs
pop <- creating.trait(pop, n.additive = 1000, 
                      trait.name = "AdditiveLargeTrait",
                      mean.target = 100, var.target = 5)

# Adding traits with 500 purely additive QTLs and 500 dominant QTLs of equal size
pop <- creating.trait(pop, n.equal.additive = 500, 
                      n.equal.dominant = 500, 
                      trait.name = "AdditiveDominantTrait",
                      mean.target = 100, var.target = 5)

# checking effects of QTLs
get.qtl.effects(pop)

# alternative to later standardize the trait
pop <- bv.standardization(pop, mean.target = 100, var.target = 5, gen=1)

pop_store <- pop
pop <- pop_store

# check TBVs
get.bv(pop, gen = 1)


# check number of generations
length(pop$breeding)

get.cohorts(pop)
# Generate 100 offspring by random mating between individuals from the two gene pools
# NOTE: Here, I do not do RANDOM mating. Random is only between the pools, i.e.,
# a random animal is drawn from pool1 as the father and a random animal is drawn
# from pool2 as mother

# PLANT BREEDERS (they don't care about sex)
pop <- breeding.diploid(pop, 
                        breeding.size = 100, 
                        selection.m.cohorts = paste0("Pool1_Founder", c("_M", "_F")), # animals from Pool1 will be used as paternal parent
                        selection.f.cohorts = paste0("Pool2_Founder", c("_M", "_F")) # animals from Pool2 will be used as maternal parent
                        )

# mating between pools is at random which means that some combinations do not occur by chance 
# to have all possible combinations, use breeding.all.combination = TRUE
# to have an individual contribute at most 20 times, set max.offspring = 10
pop <- breeding.diploid(pop, 
                        breeding.size = 100, 
                        breeding.all.combination = TRUE, 
                        max.offspring = 10,
                        selection.m.cohorts = paste0("Pool1_Founder", c("_M", "_F")), # animals from Pool1 will be used as paternal parent
                        selection.f.cohorts = paste0("Pool2_Founder", c("_M", "_F")) # animals from Pool2 will be used as maternal parent
)


# ANIMAL BREEDERS (they DO care about sex)
pop <- pop_store
pop <- breeding.diploid(pop, 
                        breeding.size = 50, 
                        name.cohort = "F1_Pool1_M_Pool2_F",
                        selection.m.cohorts = "Pool1_Founder_M", # males from Pool1 will be used as paternal parent
                        selection.f.cohorts = "Pool2_Founder_F" # females from Pool2 will be used as maternal parent
)

pop <- breeding.diploid(pop, 
                        breeding.size = 50, 
                        name.cohort = "F1_Pool1_F_Pool2_M",
                        selection.m.cohorts = "Pool1_Founder_F", # males from Pool1 will be used as paternal parent
                        selection.f.cohorts = "Pool2_Founder_M" # females from Pool2 will be used as maternal parent
                        , add.gen = 2 # this adds the offspring from this mating 
                        )
# add.gen = 2 adds the offspring from this mating to the generation 2

# ALTERNATIVE FOR ANIMAL BREEDERS
# Here, a mating list is provided to breeding diploid
#
# Example for a mate.plan:
# Generate 3 new individuals:
# First: mate individual from generation 1, sex 1, nr 1 with individual from generation 1, sex 2, nr 1
# Second: mate individual from generation 1, sex 1, nr 1 with individual from generation 1, sex 2, nr 2
# Third: mate individual from generation 1, sex 1, nr 4 with individual from generation 1, sex 2, nr 6
mate.plan <- matrix(c(1,1,1,1,2,1,
                  1,1,1,1,2,2,
                  1,1,4,1,2,6), byrow=TRUE, ncol=6)
pop <- pop_store

id.pool1.m <- get.id(pop, cohorts = "Pool1_Founder_M", use.id = T)
id.pool1.f <- get.id(pop, cohorts = "Pool1_Founder_F", use.id = T)
id.pool2.m <- get.id(pop, cohorts = "Pool2_Founder_M", use.id = T)
id.pool2.f <- get.id(pop, cohorts = "Pool2_Founder_F", use.id = T)

mate.plan <- c()
for(i in 1:50){
  mate.plan <- rbind(mate.plan, cbind(
    get.database(pop, id = sample(id.pool1.m,1))[,-3, drop=F], 
    get.database(pop, id = sample(id.pool2.f,1))[,-3, drop=F], 
    sample(c(1,0), 1) # this is the sex
  )
  )
}

for(i in 51:100){
  mate.plan <- rbind(mate.plan, cbind(
    get.database(pop, id = sample(id.pool2.m,1))[,-3, drop=F], 
    get.database(pop, id = sample(id.pool1.f,1))[,-3, drop=F], 
    sample(c(1,0), 1) # this is the sex
  )
  )
}

mate.plan

pop <- breeding.diploid(pop, fixed.breeding = mate.plan, 
                        name.cohort = "F1")


length(pop$breeding)
get.cohorts(pop)

bv.pool1 <- get.bv(pop, cohorts = paste0("Pool1_Founder", c("_M", "_F")))
bv.pool2 <- get.bv(pop, cohorts = paste0("Pool2_Founder", c("_M", "_F")))
bv.offspring <- get.bv(pop, gen = 2)
bv.parents <- get.bv(pop, gen = 1)

# Compare the genomic values of the parents and offspring
(trait.difference <- rowMeans(bv.offspring) - rowMeans(bv.parents))
# only for trait 3, which has dominant QTLs, the offspring are different (better) than the parents
# the difference between the trait with dominance is about 12


# How often was each individual used for reproduction?
ped <- get.pedigree(pop, gen = 2)
occurence <- c(ped[,2], ped[,3])
number_of_service <- table(occurence)
hist(number_of_service)
table(number_of_service)

# For trait 3 we can we much higher genomic values of the hybrid lines
# Trait architecture with dominant effects
par(mfrow=c(3,3))
hist(bv.pool1[1,], xlim=c(94,115), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv.pool1[2,], xlim=c(94,115), main="Trait 2", xlab="genomic value")
hist(bv.pool1[3,], xlim=c(94,115), main="Trait 3", xlab="genomic value")
hist(bv.pool2[1,], xlim=c(94,115), main="", ylab="Pool 2", xlab="genomic value")
hist(bv.pool2[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv.pool2[3,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv.offspring[1,], xlim=c(94,115), main="", ylab="Hybrids", xlab="genomic value")
hist(bv.offspring[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv.offspring[3,], xlim=c(94,115), main="", xlab="genomic value")
rowMeans(bv.pool1)
rowMeans(bv.pool2)
rowMeans(bv.offspring)

# Generate a PCA for all simulated individuals
get.pca(pop, gen = 1:2)
get.pca(pop, gen = 1)
get.pca(pop, gen = 2)

