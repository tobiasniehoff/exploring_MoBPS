library(MoBPS)

pop <- creating.diploid(nindi = 100, 
                        sex.s = c(rep(1,10),rep(2,90)), 
                        chr.nr = 5, 
                        chromosome.length = rep(3,5), 
                        share.genotyped = 1, 
                        nsnp = 25000
)

pop <- creating.trait(pop, n.additive = 1000, 
                      trait.name = "AdditiveTrait", 
                      mean.target = 100
)

pop <- breeding.diploid(pop, 
                        breeding.size = c(45,45))


# here, I do phenotyping
# only female offspring should be phenotyped
generation <- 2
sex <- 2
pheno.database <- cbind(generation, sex)
pop <- breeding.diploid(pop, 
                        phenotyping.database = pheno.database,
                        heritability = 0.3)

# check if only female offspring has phenotypes
get.pheno(pop, gen = 2)

# calculate average offspring phenotype first
# if you don't provide the offspring.gen, it will look up all offspring fromall generations
# in the population list
# always, only direct offspring are considered
pop <- breeding.diploid(pop, offpheno.offspring.gen = 2, offpheno.parents.gen = 1)

get.pheno.off(pop, gen = 1)

# now do a bve
pop <- breeding.diploid(pop, bve=T, bve.input.phenotype = "off",
                        relationship.matrix = "pedigree", 
                        bve.gen = 1:2)

male.founder.db <- cbind(1,1)
male.founder.bv <- get.bv(pop, database = male.founder.db)[1,]
male.founder.bve <- get.bve(pop, database = male.founder.db)[1,]
male.founder.offpheno <- get.pheno.off(pop, database = male.founder.db)[1,]

cor(male.founder.bv, male.founder.bve)
cor(male.founder.bv, male.founder.offpheno)

plot(male.founder.bv, male.founder.bve)
plot(male.founder.bv, male.founder.offpheno)

