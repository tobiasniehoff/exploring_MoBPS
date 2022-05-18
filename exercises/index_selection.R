# MoBPS 1.8.16
library(MoBPS)

pop <- creating.diploid(chr.nr = 28, nsnp = 5000, nindi = 300)


# three traits
# milk volume
# protein kg
# fat kg

cor_matrix <- matrix(NA, nrow=3, ncol=3)
diag(cor_matrix) <- 1
cor_matrix[2,1] <- 0.5
cor_matrix[3,1] <- 0.3
cor_matrix[3,2] <- -0.4

cor_matrix[upper.tri(cor_matrix)] <- t(cor_matrix)[upper.tri(cor_matrix)]

# using some made up genetic variances
varG <- c(1000,50,300)

# using some made up heritabilities
h2 <- c(0.3,0.4,0.15)

# these are the resulting error variances
varE <- (varG/h2)-varG

# these should be the mean TBVs in the first generation
target.tbv <- c(8000,280,350)

# now we add the traits
# 1000 QTL for each trait
pop <- creating.trait(pop, n.additive = rep(1000,3), 
                      trait.name = c("Milk_kg", "Protein_kg", "Fat_kg"), 
                      mean.target = target.tbv, 
                      var.target = varG,
                      shuffle.traits = 1:3, 
                      shuffle.cor = cor_matrix)


cor(t(get.bv(pop, gen=length(pop$breeding))))
cor_matrix
# the correlation is not perfectly what it should be. Increase animal numbers to get values
# closer to what is desired

# check if mean TBV is equal to target 
rowMeans(get.bv(pop, gen = length(pop$breeding)))
target.tbv

# check if variance of TBV is equal to target variances
apply(get.bv(pop, gen = length(pop$breeding)), MARGIN = 1, var)
varG

# visualize dependency of protein kg and milk kg
plot(get.bv(pop, gen = length(pop$breeding))[1,],get.bv(pop, gen = length(pop$breeding))[2,],
     ylab = "Milk_kg", xlab = "Protein_kg")
abline(lm(get.bv(pop, gen = length(pop$breeding))[2,] ~ get.bv(pop, gen = length(pop$breeding))[1,]))

# build up some LD
for(i in 1:10){
  pop <- breeding.diploid(pop, selection.size = c(20,60), breeding.size = 300)
}

#get.effective.size(pop, gen = 1)
#get.effective.size(pop, gen = length(pop$breeding))

pop_store <- pop
pop <- pop_store

# phenotyping only females
pop <- breeding.diploid(pop, 
                        phenotyping.cohorts = paste0("Cohort_",length(pop$breeding),"_F"), 
                        sigma.e = sqrt(varE))

real.varE <- apply(get.pheno(pop, cohorts = paste0("Cohort_",length(pop$breeding),"_F")),MARGIN=1,var)
real.varG <- apply(get.bv(pop, cohorts = paste0("Cohort_",length(pop$breeding),"_F")),MARGIN=1,var)

# this is the actual heritability
real.varG/(real.varG+real.varE)

# now do a PBLUP
pop <- breeding.diploid(pop, 
                        bve.gen = length(pop$breeding), 
                        depth.pedigree = 7, 
                        bve=T, 
                        relationship.matrix = "pedigree", 
                        forecast.sigma.g = F, 
                        sigma.g = sqrt(varG), 
                        sigma.e = sqrt(varE))
# Here, the variances to use for bve are directly specified to MoBPS. By default,
# MoBPS would use the variances based on TBVs as a proxy. While that is the best 
# solution, it is not what would be done in practice because variances are not 
# reestimated every generation. Here, I assume that varG and varE are estimated once
# and are reused for every bve. This is a valid assumption as for a
# reasonably big breeding program teh variance probably won't chance much over generations.
# Also, as long as the avriances are just somewhat correct, bve is not affected much.
# So this is nothin to worry about.

# specifying weights for selection index
sel.index <- c(0.5,0.3,0.2)

ebvs <- get.bve(pop, gen = length(pop$breeding), use.id = T)

# for index selection, typically the ebvs are expressed in units of 
# standard deviation, more precisely, standard deviation calculated based on
# ebvs
ebvs <- (ebvs - rowMeans(ebvs))/sqrt(diag(var(t(ebvs))))

index_values <- colSums(sel.index*ebvs, na.rm = T)

# the best 20 males should be selected
id.m <- unname(get.id(pop, cohorts = paste0("Cohort_",length(pop$breeding),"_M")))
index_values.m <- index_values[names(index_values) %in% id.m]
id.sel.m <- as.integer(names(index_values.m[order(index_values.m, decreasing = T)[1:20]]))
db.sel.m <- get.database(pop, id = id.sel.m)

# the best 60 females should be selected
id.f <- unname(get.id(pop, cohorts = paste0("Cohort_",length(pop$breeding),"_F")))
index_values.f <- index_values[names(index_values) %in% id.f]
id.sel.f <- as.integer(names(index_values.f[order(index_values.f, decreasing = T)[1:60]]))
db.sel.f <- get.database(pop, id = id.sel.f)

# breeding a new generation based on selected parents
pop.manual <- breeding.diploid(pop, 
                            breeding.size = 300, 
                            selection.m.database = db.sel.m, 
                            selection.f.database = db.sel.f)

pop.mobps <- breeding.diploid(pop, 
                         breeding.size = 300, 
                         selection.size = c(20,60),
                         multiple.bve.weights.f = sel.index, 
                         multiple.bve.weights.m = sel.index, 
                         multiple.bve.scale.m = "bve", # important
                         multiple.bve.scale.f = "bve",# important
                         selection.criteria = c("bve", "bve"),
                         selection.f.cohorts = paste0("Cohort_",length(pop$breeding),"_F"),
                         selection.m.cohorts = paste0("Cohort_",length(pop$breeding),"_M"))

ped <- get.pedigree(pop.mobps, gen = length(pop.mobps$breeding), id = T)

# if this is TRUE, the manual selection is identical to MoBPS selection
sum((sort(unique(ped[,2])) %in% sort(id.sel.m))) == length(sort(id.sel.m))

# if you would like to transform the ebvs differently, you can do that 
# manually. Selection can then either be done manually or by MoBPS
# if you insert you modified value with this function
# insert.bve()

