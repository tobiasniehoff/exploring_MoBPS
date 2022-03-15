# Mobps requires some packages
# though it could run without them, it is faster when they are used
# Torsten sent me the packages 
# it is important to have the right version of the package for the right version of 
# Mobps as miraculix is causing trouble otherwise

# for Windows
#install.packages("RandomFieldsUtils_0.6.6.tar.gz", repos = NULL, type = "source")
#install.packages("miraculix_1.0.0.1.tar.gz", repos = NULL, type = "source")

# for Linux
#install.packages("RandomFieldsUtils_1.0.6.tar.gz", repos = NULL, type = "source")
#install.packages("miraculix_1.0.5.tar.gz", repos = NULL, type = "source")

#install.packages("MoBPS_1.8.01.tar.gz", repos = NULL, type = "source")

library(MoBPS)
library(stringr)

prop_usage_external_semen <- 0.1 
# this is the proportion of semen coming from 
# other population. I will implement it so that always this proportion is 
# from the other population. In reality, you might want to say that up to 25%
# are from the other population if they have superior animals. 
# But that is a bit more tricky to code.
initial_h2 <- 0.3
n_sel_males <- 20
n_sel_females <- 25
population_size <- 200

pop <- creating.diploid(nsnp = 1000, # number of SNPs
                        chr.nr = 30, # number of chromosomes
                        nindi = population_size, # population size
                        n.additive = 100, # number of additive QTL 
                        base.bv = 100, # telling that the average true breeding value should be 100
                        share.genotyped = 0 # no animal is genotyped
                        )

# var_e <- (var_a/h2)-var_a
(var_a <- var(get.bv(pop, gen = 1)[1,]))
(var_e <- (var_a/initial_h2)-var_a)

# build up some LD by random mating with complete replacement of parent generation
# I want to build up some LD so that breeding value estimation later makes sense,
# i.e., there is some correlation/linkage between alleles and QTLs
for(i in 1:4){pop <- breeding.diploid(pop, breeding.size = population_size)}

# some generations in which the population is improved
for(i in 1:15){
  # Here, phenotyping is performed on females only
  pop <- breeding.diploid(pop,
                          sigma.e = sqrt(var_e), 
                          phenotyping.cohorts = get.cohorts(pop)[c(F,T)][length(pop$breeding)]
  )
  
  # Use at most the last 4 generations in bve. This is an arbitrary number chosen for no particular reason.
  bve.gen <- if(length(pop$breeding) >= 4){(length(pop$breeding)-3):length(pop$breeding)}else{1:length(pop$breeding)}
  
  # Here, ssGBLUP is peformed
  pop <- breeding.diploid(pop, bve = TRUE, singlestep.active = T,
                          sequenceZ = T,
                          bve.solve = "pcg",
                          bve.gen = bve.gen, 
                          sigma.e = sqrt(var_e))
  
  # Here, a new generation is created
  # Here, 10% (random) are genotyped. Note that the array used by default 
  # for genotyping includes all SNPs, so also SNPs of real QTLs
  pop <- breeding.diploid(pop, 
                          share.genotyped = 0.1, # 10% of the offspring will be genotyped randomly
                          # selection of male parents is based on bve, for females based on phenotype
                          selection.criteria = c('bve', 'pheno'), 
                          selection.size = c(n_sel_males, n_sel_females), 
                          # this defines how many males and females are to be created in the offspring generation
                          breeding.size = c(population_size/2,population_size/2),
                          avoid.mating.fullsib = TRUE,
                          avoid.mating.halfsib = TRUE,
                          selection.m.cohorts = get.cohorts(pop)[c(T,F)][length(pop$breeding)],
                          selection.f.cohorts = get.cohorts(pop)[c(F,T)][length(pop$breeding)],
                          delete.same.origin = T)
}

pop_store <- pop
pop <- pop_store

length(pop$breeding)
get.cohorts(pop)
last_common_gen <- max(as.integer(str_extract(str_extract(get.cohorts(pop), pattern = "Cohort_[[:digit:]]+"), "[[:digit:]]+")))

####
herdA.bv <- c()
herdB.bv <- c()
Fst.A.store <- c()
Fst.B.store <- c()

# Here, two populations are formed.
# 100 offspring are randomly assigned to herdA and 100 to herdB
pop <- breeding.diploid(pop,
                        share.genotyped = 0.1, # 10% of the offspring will be genotyped randomly
                        # now also the females are selected based on bve and not on phenotype anymore (for no particular reason)
                        selection.criteria = c('bve', 'bve'), 
                        selection.size = c(n_sel_males, n_sel_females),
                        breeding.size = c(population_size/2,population_size/2), 
                        avoid.mating.fullsib = TRUE,
                        avoid.mating.halfsib = TRUE,
                        name.cohort = paste0("HerdA_",last_common_gen),
                        selection.m.cohorts = get.cohorts(pop)[c(T,F)][length(pop$breeding)],
                        selection.f.cohorts = get.cohorts(pop)[c(F,T)][length(pop$breeding)],
                        delete.same.origin = T)

pop <- breeding.diploid(pop,
                        share.genotyped = 0.1, # 10% of the offspring will be genotyped randomly
                        selection.criteria = c('bve', 'bve'),
                        selection.size = c(n_sel_males, n_sel_females),
                        breeding.size = c(population_size/2,population_size/2),
                        avoid.mating.fullsib = TRUE,
                        avoid.mating.halfsib = TRUE,
                        name.cohort = paste0("HerdB_",last_common_gen),
                        selection.m.cohorts = get.cohorts(pop)[c(T,F)][length(pop$breeding)],
                        selection.f.cohorts = get.cohorts(pop)[c(F,T)][length(pop$breeding)],
                        delete.same.origin = T)

herdB.bv <- c(herdB.bv, mean(get.bv(pop, cohorts = rev(get.cohorts(pop))[1:2])))
herdA.bv <- c(herdA.bv, mean(get.bv(pop, cohorts = rev(get.cohorts(pop))[3:4])))

#######
# here, Fst values are calculated
# the definition I use is according to Wright's F-Statistic
# Fst: Expected inbreeding in subpopulations (s) 
# relative to the total population (T) they are part of
geno.B <- get.geno(pop, cohorts = rev(get.cohorts(pop))[1:2])
geno.A <- get.geno(pop, cohorts = rev(get.cohorts(pop))[3:4])

geno.B[1:5,1:6]
AF.B <- rowMeans(geno.B)/2
AF.A <- rowMeans(geno.A)/2

AF.tot <- ((AF.B+AF.A)/2)
He.totpop <- 2*AF.tot*(1-AF.tot) # expected heterozygousity in total population
He.subpop.B <- 2*AF.B*(1-AF.B) # expected heterozygousity in herd A
He.subpop.A <- 2*AF.A*(1-AF.A) # expected heterozygousity in herd B

Fst.B <- mean((He.totpop - He.subpop.B)/He.totpop, na.rm = T)
Fst.A <- mean((He.totpop - He.subpop.A)/He.totpop, na.rm = T)

Fst.A.store <- c(Fst.A.store, Fst.A)
Fst.B.store <- c(Fst.B.store, Fst.B)
#######

# 15 generations of selection
for(i in 1:15){
  cohort_names <- get.cohorts(pop)
  
  # Here, it is phenotyped
  pop <- breeding.diploid(pop,
                          sigma.e = sqrt(var_e), 
                          phenotyping.cohorts = rev(cohort_names[c(F,T)])[1:2]
  )

  # use at most the last 4 generations in bve
  # Note: The way MoBPS uses generation is not the same as we think of generations
  current_generation <- max(as.integer(str_extract(cohort_names, "[[:digit:]]+")))
  pos <- as.integer(
    str_extract_all(cohort_names, 
                    pattern = '[[:digit:]]+', 
                    simplify = T)[,1]) %in% (current_generation-3):current_generation
  
  bve.database <- get.database(pop, cohorts = cohort_names[pos])
  
  # Here, ssGBLUP is peformed
  pop <- breeding.diploid(pop, bve = TRUE, singlestep.active = T,
                          sequenceZ = T,
                          bve.solve = "pcg",
                          bve.database = bve.database,
                          sigma.e = sqrt(var_e))
  
  n_males_sel_focal_herd <- round((1-prop_usage_external_semen)*n_sel_males)
  n_males_sel_external_herd <- n_sel_males - n_males_sel_focal_herd
  
  herdA.m.bve <- get.bve(pop, cohorts = rev(cohort_names[c(T,F)])[2], use.id = T)
  herdB.m.bve <- get.bve(pop, cohorts = rev(cohort_names[c(T,F)])[1], use.id = T)
  
  new_sires_herdA <- c(
    colnames(herdA.m.bve)[order(herdA.m.bve, decreasing = T)[1:n_males_sel_focal_herd]],
    colnames(herdB.m.bve)[order(herdB.m.bve, decreasing = T)[1:n_males_sel_external_herd]]
    )
  
  new_sires_herdB <- c(
    colnames(herdB.m.bve)[order(herdB.m.bve, decreasing = T)[1:n_males_sel_focal_herd]],
    colnames(herdA.m.bve)[order(herdA.m.bve, decreasing = T)[1:n_males_sel_external_herd]]
  )
  
  # generating a new generation in herd A
  pop <- breeding.diploid(pop,
                          share.genotyped = 0.1, # 10% of the offspring will be genotyped randomly
                          selection.criteria = c('bve', 'bve'), # selection is based on bve
                          selection.size = c(n_sel_males, n_sel_females),
                          breeding.size = c(population_size/2,population_size/2),
                          avoid.mating.fullsib = TRUE,
                          avoid.mating.halfsib = TRUE,
                          name.cohort = paste0("HerdA_",current_generation+1),
                          selection.m.database = get.database(pop, id = new_sires_herdA),
                          selection.f.cohorts = rev(get.cohorts(pop)[c(F,T)])[2],
                          delete.same.origin = T)
  
  # generating a new generation in herd A
  pop <- breeding.diploid(pop,
                          share.genotyped = 0.1, # 10% of the offspring will be genotyped randomly
                          selection.criteria = c('bve', 'bve'), # selection is based on bve
                          selection.size = c(n_sel_males, n_sel_females),
                          breeding.size = c(population_size/2,population_size/2),
                          avoid.mating.fullsib = TRUE,
                          avoid.mating.halfsib = TRUE,
                          name.cohort = paste0("HerdB_",current_generation+1),
                          selection.m.database = get.database(pop, id = new_sires_herdB),
                          selection.f.cohorts = rev(get.cohorts(pop)[c(F,T)])[1],
                          delete.same.origin = T)
  
  herdB.bv <- c(herdB.bv, mean(get.bv(pop, cohorts = rev(get.cohorts(pop))[1:2])))
  herdA.bv <- c(herdA.bv, mean(get.bv(pop, cohorts = rev(get.cohorts(pop))[3:4])))
  
  #######
  # here, Fst values are calculated
  # the definition I use is according to Wright's F-Statistic
  # Fst: Expected inbreeding in subpopulations (s) 
  # relative to the total population (T) they are part of
  geno.B <- get.geno(pop, cohorts = rev(get.cohorts(pop))[1:2])
  geno.A <- get.geno(pop, cohorts = rev(get.cohorts(pop))[3:4])
  
  AF.B <- rowMeans(geno.B)/2
  AF.A <- rowMeans(geno.A)/2
  
  AF.tot <- ((AF.B+AF.A)/2)
  # next, i do marker filtering to only keep marker with at least X% minor allele 
  # frequency
  MAF_threshold <- 0.05
  AF.A <- AF.A[AF.tot > MAF_threshold & AF.tot < (1-MAF_threshold)]
  AF.B <- AF.B[AF.tot > MAF_threshold & AF.tot < (1-MAF_threshold)]
  AF.tot <- AF.tot[AF.tot > MAF_threshold & AF.tot < (1-MAF_threshold)]
  
  He.totpop <- 2*AF.tot*(1-AF.tot) # expected heterozygousity in total population (2*p*q)
  He.subpop.B <- 2*AF.B*(1-AF.B) # expected heterozygousity in herd A
  He.subpop.A <- 2*AF.A*(1-AF.A) # expected heterozygousity in herd B
  
  Fst.B <- mean((He.totpop - He.subpop.B)/He.totpop, na.rm = T)
  Fst.A <- mean((He.totpop - He.subpop.A)/He.totpop, na.rm = T)
  
  Fst.A.store <- c(Fst.A.store, Fst.A)
  Fst.B.store <- c(Fst.B.store, Fst.B)
  #######
}

# plotting trajectory if mean true breeding values
plot(herdA.bv, col = "dodgerblue", pch = 20, ylab = "TBV", xlab = "generation")
points(herdB.bv, col = "tomato", pch = 20)
legend(legend=c("Average tbv herd A", "Average tbv herd B"), 
      fill = c('dodgerblue', 'tomato'), "topleft")


# plotting trajectory if mean true breeding values
plot(Fst.A.store, col = "dodgerblue", pch = 20, ylab = "TBV", xlab = "generation")
points(Fst.B.store, col = "tomato", pch = 20)
legend(legend=c("Fst herd A", "Fst herd B"), 
       fill = c('dodgerblue', 'tomato'), "topleft")
