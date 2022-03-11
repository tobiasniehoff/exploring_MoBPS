# for Windows
#install.packages("MoBPS_1.8.01.tar.gz", repos = NULL, type = "source")
#install.packages("PRandomFieldsUtils_0.6.6.tar.gz", repos = NULL, type = "source")
#install.packages("miraculix_1.0.0.1.tar.gz", repos = NULL, type = "source")

# for Linux
#install.packages("MoBPS_1.8.01.tar.gz", repos = NULL, type = "source")
#install.packages("RandomFieldsUtils_1.0.6.tar.gz", repos = NULL, type = "source")
#install.packages("miraculix_1.0.5.tar.gz", repos = NULL, type = "source")

library(MoBPS)

# simulate historic population as in Jibrila et al. 2020
n.snps <- 60000
n.chr <- 30
n.qtls <- 3000
# because later I don't want real QTLs to be on SNPs that are on the array
n.snps <- n.snps+n.qtls 
# drawing effects from a gamma distribution with shape parameter 0.4
effect0 <- effect2 <- qtl.effects <- rgamma(n = n.qtls, shape = 0.4)

qtl.snps <- seq(1, (n.snps), (n.snps)/(n.qtls)) # these are the SNPs where I will put a QTL
# the positions are evenly distributed

SNP <- qtl.snps
SNP <- (SNP-1)%%(n.snps/n.chr)+1

pos <- sample(x = n.qtls, size = n.qtls/2,replace = F)
effect0[pos] <- qtl.effects[pos]
effect2[pos] <- 0
effect0[-pos] <- 0
effect2[-pos] <- qtl.effects[-pos]
effect1 <- rowMeans(cbind(effect0, effect2))
chromosome <- sort(rep(1:n.chr,n.qtls/n.chr))
effect.matrix <- cbind(SNP, chromosome, effect0, effect1, effect2)

genotyping_array <- rep(F, n.snps)
genotyping_array[-qtl.snps] <- T

# actually the marker array and effect matrix stuff above does not matter 
# for the simulation of pop history
# as no genotyping or selection based on performance is conducted
################################################################################

ninds <- 5000
population <- creating.diploid(nsnp = n.snps, nindi = ninds, "random",
                               chr.nr = n.chr, chromosome.length = 1,
                               real.bv.add = effect.matrix, snps.equidistant = T, 
                               name.cohort = "Founder")

population$info$array.name # see array names
population <- add.array(population,marker.included = genotyping_array)
# This adds the marker array. This array is the 2nd array as the first one
# always includes all SNPs
population$info$array.name # see array names

# simulating 2997 generations of random mating with linearly decreasing
# population size from 5000 to 50
d <- data.frame(x=c(1,2997),y=c(5000,50))
l <- lm(y~x, d)

gen <- c()
for(i in 1:2997){
  ninds <- round(predict(l, rbind(d,c(i,NA)))[[3]])
  ninds <- ninds + (ninds%%2) # this line makes sure that I only have even
  # numbers so that the sex ratio is exactly 1:1
  population <- breeding.diploid(population, breeding.size = ninds)
  
  # defining a new base generation and deleting old data is more memory efficient 
  # and thus faster on my machine
  # This has an effect after some high 10s or 100s generations.
  population <- new.base.generation(population, 
                                    base.gen = length(population$breeding), 
                                    delete.previous.gen = T, 
                                    delete.bve.data = T, delete.breeding.totals = T)
  gen <- i
}
gen

# population expansion to 5000 inds in 3 generations
d <- data.frame(x=c(2997,3000),y=c(50,5000))
l <- lm(y~x, d)
gen <- c()
for(i in 2998:2999){
  ninds <- round(predict(l, rbind(d,c(i,NA)))[[3]])
  # this line makes sure that the population has a size that allows
  # for the same number of males as females (sex ratio is exactly 1:1)
  ninds <- ninds + (ninds%%2) 
  population <- breeding.diploid(population, breeding.size = ninds)
  population <- new.base.generation(population,
                                    base.gen = length(population$breeding),
                                    delete.previous.gen = T,
                                    delete.bve.data = T, delete.breeding.totals = T)
  gen <- i
}
gen

# this is the 3000th generation
population <- breeding.diploid(population, breeding.size = c(100,1000)) # these will be my founders

get.vcf(population, path = "founder_population", gen = length(population$breeding))
save.image(file = "founder_population.RData")

# With this line, you can create the populatin object based on a vcf file.
# Note that this does not load the marker effects.
population <- creating.diploid(vcf="population.vcf", bpcm.conversion = 1000000, 
                               name.cohort = "Historic")

# With these lines, you can add the marker array and effect matrix.
population <- creating.trait(population, real.bv.add = effect.matrix)
marker.array <- 1:n.snps %in% qtl.snps
population <- add.array(population, marker.included = !marker.array)

