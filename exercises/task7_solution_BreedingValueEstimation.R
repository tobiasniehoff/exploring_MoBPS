# R version 4.0.2 (2020-06-22)
# MoBPS_1.8.07

library(MoBPS)

pop <- creating.diploid(nindi = 100, 
                        sex.s = c(rep(1,50),rep(2,50)), 
                        chr.nr = 5, 
                        chromosome.length = rep(3,5), 
                        share.genotyped = 0, 
                        nsnp = 25000
                        )

# you could also specify how many marker should be on every chromosome
pop <- creating.diploid(nindi = 100, 
                        sex.s = c(rep(1,82),rep(2,18)), # this vector tells which animal should be male or female 
                        #chr.nr = 5, 
                        chromosome.length = c(5,4,3,1,1), 
                        share.genotyped = 0, 
                        nsnp = c(2000,1000,1000,500,10000,10500)
)

# get.map() is not working in MoBPS_1.8.07
map <- get.map(pop)
dim(map)

# plot the position of the markers in centi Morgan
plot(as.numeric(map[,3]))

pop <- creating.trait(pop, n.additive = 1000, 
                      trait.name = "AdditiveTrait", 
                      mean.target = 100
                      )
mean(get.bv(pop, gen = 1)[1,])


# I don't want to specify heritability later as the heritability
# is a feature of the population AND the trait. If the population changes
# due to breeding, heritability is not suitable. Better use sigma.e
# because the environmental variance does not change
target.heritability.founders <- 0.3
var_a <- var(get.bv(pop, gen = 1)[1,])
# var_e <- (var_a/h2)-var_a
var_e <- (var_a/target.heritability.founders)-var_a

# the population list should have 12 generations
# the first generation is the founder generation
# the generations are distinct. The last generation is used as parents 
# for the next one
for(i in 2:12){
  
  if(i <=10){
    # generation of new animals
    pop <- breeding.diploid(pop, 
                            breeding.size = c(50,50), 
                            share.genotyped = 0)
  }
  
  if(i >= 10){
    # only males should be phenotyped in generation 10
    if(i == 10){
      pheno.db <- get.database(pop, gen = i)
      genotyping.db <- pheno.db
      pheno.db <- pheno.db[pheno.db[,2]==1,, drop=F]
    }
    if(i == 11){
      pheno.db <- get.database(pop, gen = i)
      genotyping.db <- pheno.db
    }
    # in generation 12 all are phenotyped but only 20% are genotyped
    if(i == 12){
      pheno.db <- get.database(pop, gen = i)
      id.12 <- get.id(pop, gen = 12, use.id = T)
      genotyping.db <- get.database(pop, 
                                    id = sample(id.12, size = ceiling(length(id.12)*0.2)))
    }
    
    # here, genotyping is done
    pop <- breeding.diploid(pop, 
                            genotyped.database = genotyping.db, 
                            share.genotyped = 1)
    
    # here, phenotyping is done
    pop <- breeding.diploid(pop, 
                            phenotyping.database = pheno.db,
                            sigma.e = sqrt(var_e))
    
    # here, bve is performed
    # I use ssGBLUP
    pop <- breeding.diploid(pop, bve = TRUE,
                            singlestep.active = T, 
                            bve.gen = 10:i, 
                            sigma.e = sqrt(var_e))
    
    # the 12th generation is only evaluated but they don't have offspring
    if(i != 12){
      # here, a new generation is created
      # automatically, the previous generation is used as parents
      pop <- breeding.diploid(pop, 
                              breeding.size = c(50,50),
                              selection.criteria = c("bve", "bve"), 
                              share.genotyped = 0
      )
    }
  }
}

# Alternative solution (and slightly different breeding program):
# Here, heritability is used but sigma.e could have been provided just as well
# this only goes to show what different options are available in MoBPS
# population <- pop
# for(index in 1:10){
#   population <- breeding.diploid(population, breeding.size = 100)
# }
# # In generation 10 only males are phenotyped
# population <- breeding.diploid(population, phenotyping.database = cbind(10,1), heritability = 0.3)
# 
# # On default all individuals are assumed to be genotyped
# # use share.genotyped to generate some non-genotyped individuals
# population <- breeding.diploid(population, breeding.size = 100, share.genotyped = 0.2)
# 
# # Generation of phenotypic data
# # In generation 11 & 12 all individuals are phenotyped
# population <- breeding.diploid(population, phenotyping.gen = 11:12, heritability = 0.3)
# 
# pop <- population



# if you didn't manage to set up the program, load in the data
#load("population.RData")
#pop <- population

# checking which animal was genotyped in generation 10
get.genotyped(pop, gen = 10)
# checking which animal was genotyped in generation 12
get.genotyped(pop, gen = 12)
length(pop$breeding)
get.cohorts(pop)

pop_store <- pop
pop <- pop_store


### subpart 2
# To only base bve on a subset of markers, they are randomly sampled first
# and then an array with only the subset is created in MoBPS. The user can 
# then specify the array to use

num_snps <- nrow(get.map(pop))
# alternatively: num_snps <- nrow(get.geno(pop, gen =1))
n_random_snps <- 100
resample <- function(x, ...) x[sample.int(length(x), ...)]
included.marker <- resample(c(rep(T, n_random_snps), rep(F, num_snps-n_random_snps)))
sum(included.marker)

# Alternatively:
# included.marker <- rep(FALSE, 25000)
# included.marker[sample(1:25000, 100)] <- TRUE

pop <- add.array(pop, marker.included = included.marker, array.name = "array100")

num_snps <- nrow(get.map(pop))
n_random_snps <- 10000
included.marker <- resample(c(rep(T, n_random_snps), rep(F, num_snps-n_random_snps)))
sum(included.marker)

pop <- add.array(pop, marker.included = included.marker, array.name = "array10000")

# these arrays are defined for the population
pop$info$array.name


# single-step bve
# single-step is a combination of a pedigree matrix and a marker matrix
# which is useful when only part of the population is genotyped but
# non genotyped animals with pedigree records should still be used in bve
pop <- breeding.diploid(pop, bve = T, 
                        bve.gen = 10:12,
                        singlestep.active = T, 
                        bve.array = "array100", 
                        sigma.e = sqrt(var_e)
                        )

pop <- breeding.diploid(pop, bve = T, 
                        bve.gen = 10:12,
                        singlestep.active = T, 
                        bve.array = "array10000", 
                        sigma.e = sqrt(var_e)
)

# bve.array = 3 is identical to bve.array = "array10000"
# remove.effect.position = TRUE means that in bve, SNPs positions that have a 
# QTL effect are not used as this may be a string assumption and not reflect
# reality well
pop <- breeding.diploid(pop, bve = T, 
                        bve.gen = 10:12,
                        singlestep.active = T, 
                        bve.array = 3, 
                        remove.effect.position = TRUE,
                        sigma.e = sqrt(var_e)
)

plot(get.bv(pop, gen =10)[1,], get.bve(pop, gen = 10)[1,])

analyze.bv(pop, database = cbind(10,2))

# pedigree BLUP
pop <- breeding.diploid(pop, bve = T, 
                        bve.gen = 10:12, 
                        relationship.matrix = "pedigree", 
                        bve.array = "array100", 
                        sigma.e = sqrt(var_e)
)

# GBLUP
pop <- breeding.diploid(pop, bve = T, 
                        bve.gen = 10:12, 
                        relationship.matrix = "vanRaden", 
                        bve.array = "array100", 
                        sigma.e = sqrt(var_e)
)
