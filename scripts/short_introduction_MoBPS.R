
# MoBPS_1.11.46 
library(MoBPS)

# species with 30 chromosomes, 8000 loci, 20 founder individuals
# and two correlated traits. The correlation between the traits in 0.3.
# The first trait has 10 additive QTLs, the second one has 20
pop <- creating.diploid(nindi = 2000, chr.nr = 30, chromosome.length = 1, nsnp = 8000, 
                        n.additive = c(1000,2000)  
                        , trait.cor = matrix(c(1,0.3,0.3,1),nrow=2,ncol=2)
                        )

map <- get.map(pop)

get.id(pop, gen = 1)
get.id(pop, cohorts = 'Cohort_1_M')
db <- get.database(pop, id = c(1,2,3,5,15))

# this function extracts the true breeding values of the animals
get.bv(pop, database = db, use.id = TRUE)

geno <- get.geno(pop, database = db)

# can also use the id to get databases
get.bv(pop, gen = TRUE, use.id = TRUE)[,order(get.bv(pop, gen = TRUE, use.id = TRUE))]
# say we wanted to get the datbase of the animals with the two lowst breeding values
# which are iD 7 and ID 3 here
get.database(pop, id = c(7,3))

dim(geno)
geno[1:10, ]

# extracts the TRUE qtl effects
q <- get.qtl.effects(pop)

# these are the true additive QTL effects
# the first list is effects for trait 1 and the second list is for trait 2
q[[1]]
# in creating.diploid(), we said the first trait should have 10 QTLs and the second one
# should have 20. But we see that the second one has 30.
# this is because we what the two traits to be correlated. If we want them to be correlated
# MoBPS uses the first 10 SNPs of trait 2 to be correlated to SNPs of trait 1
# and the last 20 (or less but never less than asked for in input)
# are not correlated to trait1 SNPs but are needed for the required variance
qtl.trait1 <- q[[1]][[1]]
qtl.trait2 <- q[[1]][[2]]

# these are the true dominance effects. We haven't simulated any so there ar none
q[[2]]


# achieved correlation between breeding values
tbv.trait1 <- get.bv(pop, gen = get.ngen(pop), use.id = TRUE)[1,]
tbv.trait2 <- get.bv(pop, gen = get.ngen(pop), use.id = TRUE)[2,]
cor(tbv.trait1, tbv.trait2)
# this correlation is different from what we asked for. This is
# because of the low number of QTL and individuals which makes it
# difficult to tune the effects properly. 
# Try creating a founder population again and simulate 2000 founders and
# 1000 and 2000 QTLs and see what happens

# correlation of SNP effects at loci that affect both traits
pleiotropic.snp <- qtl.trait2[,6] %in% qtl.trait1[,6]

cor(qtl.trait1[,5], qtl.trait2[pleiotropic.snp,5])

# number of pleiotropic QTL
sum(pleiotropic.snp)
# number of QTL affecting trait 2
length(qtl.trait2[,1])

################################################################################
################################################################################
# here, we adjust genetic variance to our needs

# 37 chromosomes of random length
pop <- creating.diploid(nindi = 20, 
                        chr.nr = 30, 
                        chromosome.length = abs(runif(37,min = 0.5, max = 2)), 
                        nsnp = 8000, 
                        n.additive = 2000
)

# here, we will adjust the genetic variance to our needs
# this could have also been done before in creating.diploid()
# we will 

var.before <- var(get.bv(pop, gen = 1)[1,])
mean.before <- mean(get.bv(pop, gen = 1)[1,])

q <- get.qtl.effects(pop)
qtl.before <- q[[1]][[1]]
additive.eff.before <- (qtl.before[,5] - qtl.before[,3])/2

# this function adjusts the variance and mean of the trait
# You can also set this in "creating.diploid()"
# But sometimes it is easier to do it later, or sometimes
# you would like to re-adjust the mean, variance and correlation
pop <- bv.standardization(pop, var.target = 1,gen = 1
                          , traits = 1
                          , mean.target = 42)

var.after <- var(get.bv(pop, gen = 1)[1,])
mean.after <- mean(get.bv(pop, gen = 1)[1,])

q <- get.qtl.effects(pop)
qtl.after <- q[[1]][[1]]
additive.eff.after <- (qtl.after[,5] - qtl.after[,3])/2

# see that the genetic variance was changed
# MoBPS does this by changing the QTL effects
var.before
var.after

# average TBV change
# Note that "TBV" is MoBPS is not a deviation from an average but rather the true
# performance the animal has
mean.before
mean.after

# the QTL effects have been re-adjusted
mean(additive.eff.before)
mean(additive.eff.after)

# but they are perfectly correlated to teh effects before adjustment
# meaning that they have only been scaled
cor(additive.eff.before, additive.eff.after)
plot(additive.eff.before, additive.eff.after)

###############################
# generate new animals

# this generates new offspring
# sires are the first 3 males and all females
# mating is random
# it will create 5 male and 20 female offspring
pop <- breeding.diploid(pop, breeding.size = c(5,20), 
                        selection.m.database = cbind(1,1,1,3), 
                        selection.f.cohorts = "Cohort_1_F")

# you see that new cohorts have been created for the new animals
get.cohorts(pop)

# you see that 5 males and 20 females have been created
get.database(pop, gen = 2)

# this shows you the pedigree of the 5 newly created males
get.pedigree(pop, database = cbind(2,1,1,5))

# you can also name the new cohorts
# also,we say that we can have at most 5 offspring from a particular sirexdam combination
# EVery individual can have at most 7 offspring
pop <- breeding.diploid(pop, breeding.size = 34, name.cohort = "OUR_COHORT_name"
                        , max.mating.pair = 5, max.offspring = 7)

# by default, _M and _F is always added
get.cohorts(pop)

# this returns the allele frequency for our cohorts
get.allele.freq(pop, cohorts = c("OUR_COHORT_name_M", "OUR_COHORT_name_F"))

################################################################################
# short example for breeding value estimation

pop <- creating.diploid(nindi = 10, nsnp = 5000, chr.nr = 30, chromosome.length = 1, n.additive = 1000)

# generate some LD by random mating(though not necessarily needed for my point)
# so that genomics will work later
for(i in 1:10){
  pop <- breeding.diploid(pop, breeding.size = 200, delete.same.origin = TRUE)
}
for(i in 1:10){
  pop <- breeding.diploid(pop, breeding.size = 200, delete.same.origin = TRUE, selection.size = c(20,20))
}


varA <- var(get.bv(pop, gen = 1)[1,])
# this is the h2 we want
target.h2 <- 0.3
varE <- (varA/target.h2)- varA

varA / (varA+ varE)

# here we phenotype the animals
# with sigma.e, we provide the error (environmental) standard deviation
pop <- breeding.diploid(pop, phenotyping.gen = 1:get.ngen(pop), sigma.e = sqrt(varE))

# here we genotype the youngest animals. By default all loci will be genotyped (like sequencing)
pop <- breeding.diploid(pop, genotyped.gen = get.ngen(pop))

# here we perform pedigree breeding value estimation
pop <- breeding.diploid(pop, bve = TRUE, 
                        relationship.matrix = "pedigree", 
                        sigma.e = sqrt(varE), 
                        sigma.g = sqrt(varA), 
                        bve.gen = 19:get.ngen(pop))
ebv.PBLUP <- get.bve(pop, gen = get.ngen(pop))[1,]
ebv.PBLUP.all <- get.bve(pop, gen = 19:get.ngen(pop))[1,]

# here we do GBLUP. Note that we only use the youngest animals as only generation 11 was genotyped
pop <- breeding.diploid(pop, bve = TRUE, 
                        relationship.matrix = "vanRaden", 
                        remove.effect.position = TRUE, # this option means that we do not include the positions
                        # of real QTL in construction of our G matrix as it is a unrealistic 
                        # assumption that we genotype the true QTL positions
                        sigma.e = sqrt(varE), 
                        sigma.g = sqrt(varA), 
                        bve.gen = get.ngen(pop))
ebv.GBLUP <- get.bve(pop, gen = get.ngen(pop))[1,]

# here we do single step GBLUP
pop <- breeding.diploid(pop, bve = TRUE, 
                        relationship.matrix = "vanRaden", 
                        singlestep.active = TRUE,
                        remove.effect.position = TRUE, # this option means that we do not include the positions
                        # of real QTL in construction of our G matrix as it is a unrealistic 
                        # assumption that we genotype the true QTL positions
                        sigma.e = sqrt(varE), 
                        sigma.g = sqrt(varA), 
                        bve.gen = 19:get.ngen(pop) 
                        )
ebv.ssGBLUP <- get.bve(pop, gen = get.ngen(pop))[1,]
ebv.ssGBLUP.all <- get.bve(pop, gen = 19:get.ngen(pop))[1,]

pheno <- get.pheno(pop, gen = get.ngen(pop))[1,]
pheno.all <- get.pheno(pop, gen = 19:get.ngen(pop))[1,]
tbv <- get.bv(pop, gen = get.ngen(pop))[1,]
tbv.all <- get.bv(pop, gen = 19:get.ngen(pop))[1,]

# the accuracy for all animals in the training set
cor(tbv.all,pheno.all)
cor(tbv.all,ebv.PBLUP.all)
cor(tbv.all,ebv.ssGBLUP.all)

# check the prediction accuracies for the youngest animals
cor(tbv,pheno)
cor(tbv,ebv.PBLUP)
cor(tbv,ebv.GBLUP)
cor(tbv,ebv.ssGBLUP)


# does the same
analyze.bv(pop, gen = 21)
analyze.bv(pop, gen = 19:21)

################################################################################
pop <- creating.diploid(nindi = 20, nsnp = 25000, chr.nr = 30, chromosome.length = 1, n.additive = 5000)

# generate some LD by random mating(though not necessarily needed for my point)
# so that genomics will work later
for(i in 1:10){
  pop <- breeding.diploid(pop, breeding.size = 200, delete.same.origin = TRUE)
}

varA <- var(get.bv(pop, gen = 1)[1,])
# this is the h2 we want
target.h2 <- 0.3
varE <- (varA/target.h2)- varA

varA / (varA+ varE)

# for several generations, a simple breeding program using GBLUP could look like this
# I already use overlapping generations here
# The males can be selected from this generation and the previous generation
# females can be selected from this generation and 3 generations before

for(i in 1:10){
  pop <- breeding.diploid(pop, phenotyping.gen = get.ngen(pop), sigma.e = sqrt(varE))
  pop <- breeding.diploid(pop, genotyped.gen = get.ngen(pop))
  pop <- breeding.diploid(pop, bve = TRUE, 
                          relationship.matrix = "vanRaden", 
                          remove.effect.position = TRUE, # this option means that we do not include the positions
                          # of real QTL in construction of our G matrix as it is a unrealistic 
                          # assumption that we genotype the true QTL positions
                          sigma.e = sqrt(varE), 
                          sigma.g = sqrt(varA), 
                          bve.gen = 11:get.ngen(pop))
  
  # select 20 males and 50 females based on their EBV
  pop <- breeding.diploid(pop, selection.criteria = 'bve', 
                          selection.size = c(20,50), 
                          breeding.size = 200, 
                          delete.same.origin = TRUE
                          , selection.m.cohorts = paste0("Cohort_",(get.ngen(pop)-1):get.ngen(pop),"_M" )
                          , selection.f.cohorts = paste0("Cohort_",(get.ngen(pop)-3):get.ngen(pop),"_F" )
                          )
}
get.cohorts(pop)
# by checking the pedigree, you can see that animals from different generation were selected
get.pedigree(pop, gen = get.ngen(pop))

av.genetic.level <- c()
for(i in 11:get.ngen(pop)){
  av.genetic.level <- c(av.genetic.level, mean(get.bv(pop, gen = i)))
}
av.genetic.level <- av.genetic.level - av.genetic.level[1]


av.F.level <- c()
for(i in 11:get.ngen(pop)){
  av.F.level <- c(av.F.level, mean(inbreeding.emp(pop, gen = i)) )
}
av.F.level <- av.F.level - av.F.level[1]

mat.age.sire <- matrix(NA,nrow = 5,ncol = 11)
rownames(mat.age.sire) <- paste0(1:5, "_years")
for(i in 11:get.ngen(pop)){
  tab <- table(i - get.pedigree(pop, gen = i, raw = TRUE)[,4])
  mat.age.sire[1:length(tab),i-10] <- tab
}

mat.age.dam <- matrix(NA,nrow = 5,ncol = 11)
rownames(mat.age.dam) <- paste0(1:5, "_years")
for(i in 11:get.ngen(pop)){
  tab <- table(i - get.pedigree(pop, gen = i, raw = TRUE)[,7])
  mat.age.dam[1:length(tab),i-10] <- tab
}

plot(11:get.ngen(pop), av.genetic.level, type = "b", xlab = "generation", ylab = "Average TBV")
plot(11:get.ngen(pop), av.F.level, type = "b", xlab = "generation", ylab = "Average IBD inbreeding level")

colnames(mat.age.sire) <- 11:21
mat.age.sire <- mat.age.sire[1:2,]
mat.age.sire[is.na(mat.age.sire)] <- 0
mat.age.sire <- 100*mat.age.sire/colSums(mat.age.sire)
colosr <- c("green", "orange", "blue", "purple")

colnames(mat.age.dam) <- 11:21
mat.age.dam <- mat.age.dam[1:4,]
mat.age.dam[is.na(mat.age.dam)] <- 0
mat.age.dam <- 100*mat.age.dam/colSums(mat.age.dam)


par(mfrow = c(1, 2))

barplot(mat.age.dam,
        main = "Dams",
        xlab = "Year",
        ylab = "%",
        axes = TRUE, 
        space = 0, 
        col = colosr)

barplot(mat.age.sire,
        main = "Sires",
        xlab = "Year",
        ylab = "%",
        axes = TRUE, 
        space = 0, 
        col = colosr,
        legend.text = paste0(1:4, " years old"),
        args.legend = list(x = "bottomright"))

