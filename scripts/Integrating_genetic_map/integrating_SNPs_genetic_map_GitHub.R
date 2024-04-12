
# MoBPS_1.11.40
library(MoBPS)

# In this script, I am reading in the genetic map from Groenen et al. 2009
# https://genome.cshlp.org/content/19/3/510.short
# A high-density SNP-based linkage map of the chicken genome reveals sequence features correlated with recombination rate
#
# After my cleaning, it covers 28 autosomes. I removed their data from other linkage groups and
# sex chromosomes
# After cleaning, 8620 SNPs remain. For higher density panels, one could 
# use the physical position and extrapolate where other SNPs would be located on the genetic map
# Also, I use teh average cM position even though the authors provide sex 
# specific maps
#
# The information was in their Supp Table S1.xls (https://genome.cshlp.org/content/19/3/510/suppl/DC1)
setwd("C:/Users/nieho006/Documents/PhD Wageningen/MoBPS/chicken_genetic_map")

IN <- as.data.frame(readxl::read_xls("Supplemental_Table_S1.xls", col_names = TRUE))

head(IN)
tail(IN, n = 100)

table(IN[,1])

# subsetting only to autosomes with numbers (1 to 28)
IN <- IN[IN[,1] %in% 1:28,]

# order according to cM position
IN <- IN[order(IN[,"average (cM)"]),]

# order according to chromosome
IN <- IN[order(IN[,"Chromosome"]),]

head(IN)

# extracting the info that is relevant for MobPS
# (actually teh bp position is not relevant)
map.for.MoBPS <- IN[,c(1,2,5,6)]

# converting from centi Morgan to Morgan
map.for.MoBPS[,3] <- map.for.MoBPS[,3] / 100

table(map.for.MoBPS[,1])
dim(map.for.MoBPS)
head(map.for.MoBPS)
colnames(map.for.MoBPS) <- c("Chromosome", "SNP ID", "Morgan position", "BP position")

# MoBPS does not accept a Morgan position to be 0
# so that is why we will add a little bit on top, menaing shoft all positions by a small amount
chrom.i <- 1
for(chrom.i in 1:length(unique(map.for.MoBPS[,1]))){
  map.for.MoBPS[map.for.MoBPS[,1] == chrom.i,3] <- map.for.MoBPS[map.for.MoBPS[,1] == chrom.i,3] + 0.001
}

head(map.for.MoBPS)

any(is.na(map.for.MoBPS[,1]))
any(is.na(map.for.MoBPS[,2]))
any(is.na(map.for.MoBPS[,3]))
# there are some SNP for which only teh gentic position but not the 
# physical position is known
any(is.na(map.for.MoBPS[,4]))


pop <- creating.diploid(nindi = 100, map = map.for.MoBPS)

# the genetic map was handed over successfully
map <- get.map(pop)

len.Morgan <- c()
chrom.i <- 1
for(chrom.i in 1:length(unique(map[,1]))){
  map.chromo <- map[map[,1] == unique(map[,1])[chrom.i], ]
  len.Morgan <- c(len.Morgan, max(as.numeric(map.chromo[,3])) - min(as.numeric(map.chromo[,3])) )
}
len.Morgan
barplot(height = len.Morgan, names.arg = 1:28, xlab = "chromosome", ylab = "length in Morgan", ylim = c(0,5))

################################################################################
################################################################################

map.groenen <- get.map(pop)

# this function returns conversion functions for every chromosome
# based on an input genetic map
# there may be chromosomes in teh genome that are not on the reference genetic map
# in that case, only the bp.per.Morgan argument is used to make a conversion factor/function
# bp.per.Morgan is only used if that chromosome is not on the genetic map
#
# Update 10.04.2024: 
# It may be that the minimum and maximum physical positions of the genetic map are
# higher or lower than the positions of SNPs that we would later like to integrate. 
# Without modification, all SNPs beyond the ends of the chromosomes would get the same
# genetic positions. An alternative may be to do an linear interpolation to the actual 
# genetic length known from some other publication.
# Another option is to just use the ordinary conversion factor for those regions.
# This option is implemented in this function.
#
# Update 12.04.2024: 
# bp.per.Morgan.chrom can accept a vector of length of the chromosomes so that
# chromosome specific conversion factors can be supplied
# added argument "estimate.conversion.factor". If this is TRUE,
# and bp.per.Morgan is not a vector, conversion factors are estimated based on the 
# provided genetic map (length_bp/length_Morgan), i.e., chromosome specific conversion factors
function.bp.Morgan <- function(map, 
                               bp.per.Morgan = 30000000, 
                               bp.per.Morgan.chrom = NA,
                               n.chromosomes = 30,
                               return.cleaned.map = FALSE,
                               estimate.conversion.factor = TRUE){
  
  map <- map[order(as.numeric(map[,4])),]
  map <- map[order(as.numeric(map[,1])),]
  
  list.conversion.functions.chromosomes <- list()
  list.min.max.bp.chromosomes <- list()
  cleaned.map <- c()
  chromosome.conversion.factor <- c()
  
  chrom.i <- 1
  for(chrom.i in 1:length(unique(map[,1]))){
    pos <- map[,1] == unique(map[,1])[chrom.i]
    
    removed.snps <- 0
    
    map.chrom <- map[pos,]
    removed.snps <-  removed.snps + sum(duplicated(map.chrom[,2]))
    map.chrom <- map.chrom[!duplicated(map.chrom[,2]),]
    
    bp <- as.numeric(map.chrom[,4])
    M <- as.numeric(map.chrom[,3])
    cor(order(bp), order(M))
    
    # first, removing all SNPs that are at bp position 0
    removed.snps <- removed.snps + sum(as.numeric(map.chrom[,4]) == 0)
    map.chrom <- map.chrom[as.numeric(map.chrom[,4]) != 0,]
    
    # next, we remove SNPs that don't fit in
    repeat{
      if(nrow(map.chrom) < 10){break}
      
      bp <- as.numeric(map.chrom[,4])
      map.chrom <- map.chrom[order(bp),]
      bp <- as.numeric(map.chrom[,4])
      M <- as.numeric(map.chrom[,3])
      
      #hist(bp)
      
      rm.pos <- bp > (mean(bp) + 6*sd(bp))
      rm.pos <- rm.pos | (bp < (mean(bp) - 6*sd(bp)))
      
      rm.pos <- rm.pos | (M > (mean(M) + 6*sd(M)))
      rm.pos <- rm.pos | (M < (mean(M) - 6*sd(M)))
      
      if(!any(rm.pos)){break}
      
      map.chrom <- map.chrom[!rm.pos,]
      removed.snps <- removed.snps + sum(rm.pos)
    }
    cor(order(bp), order(M))
    
    # here, we estimate a linear regression and check which SNPs deviate the most
    # outliers will be removed
    repeat{
      if(nrow(map.chrom) < 10){break}
      
      bp <- as.numeric(map.chrom[,4])
      map.chrom <- map.chrom[order(bp),]
      bp <- as.numeric(map.chrom[,4])
      M <- as.numeric(map.chrom[,3])
      
      plot(bp, M)
      l <- lm(M ~ bp)
      abline(l)
      s <- summary(l)
      sd.res <- sd(s$residuals)
      
      ixk <- function(p) {
        x<-qnorm(1-p)		#truncation point
        i<-dnorm(x)/p		#selection intensity
        k<-i*(i-x)		#variance reduction coefficient
        return(c(i,x,k))
      }
      cor(order(bp), order(M))
      
      # decide that outliers must be at least 4 SD away from the mean
      # OR more than X standard deviations which correspond to something 
      # from selection intensity. Basically, in a very large sample, you will have some 
      # observations that deviate a lot from the mean but they are not necessarily outliers 
      dev <- (abs(s$residuals) / sd.res)
      higher.1 <- dev/max(c(4,ixk(p = 1/nrow(map.chrom))[1]))
      if(!any(higher.1 > 1)){break}
      
      rm.snp <- which.max(dev/max(c(4,ixk(p = 1/nrow(map.chrom))[1])))
      
      map.chrom <- map.chrom[-rm.snp,]
    }
    cor(order(bp), order(M))
    
    # here, I remove SNPs due to rank issues
    repeat{
      if(nrow(map.chrom) < 10){break}
      
      bp <- as.numeric(map.chrom[,4])
      map.chrom <- map.chrom[order(bp),]
      M <- as.numeric(map.chrom[,3])
      bp <- as.numeric(map.chrom[,4])
      
      plot(bp, M)
      l <- lm(M ~ bp)
      abline(l)
      summary(l)
      
      diff <- order(M, decreasing = FALSE)[-1] - order(M)[-length(M)]
      if( ( sum(diff == 1)/length(diff) ) == 1 ){break}
      
      min(which(diff != 1))
      abs(diff[diff != 1])
      
      map.chrom <- map.chrom[-(min(which(diff != 1))+1),]
      
      # an alternative that also works
      #map.chrom <- map.chrom[-which.max(abs(diff - 1)),]
      
      removed.snps <- removed.snps + 1
    }
    cor(order(bp), order(M))
    
    if(nrow(map.chrom) > 10){
      bp <- as.numeric(map.chrom[,4])
      M <- as.numeric(map.chrom[,3])
      
      list.min.max.bp.chromosomes[[chrom.i]] <- c(min(bp), max(bp))
      
      cor(order(bp), order(M))
      removed.snps <- nrow(map[pos,]) - nrow(map.chrom)
      #cat(paste0("Removed ", removed.snps, " SNPs from chromosome ", chrom.i," because of rank issues.\n"))
      cat(paste0("Removed ", removed.snps, " SNPs from chromosome ", chrom.i,".\n"))
    }
    
    if(length(unique(map.chrom[,3])) < 10 | length(unique(map.chrom[,4])) < 10){
      # I want at least 10 unique data points to estimate the function
      # this is if there is not enough data to estimate the function
      
      est.conversion.function <- NA
      
      est.conversion.function <- 1/bp.per.Morgan
      
      if(length(bp.per.Morgan.chrom) > 1){
        if(!is.na(bp.per.Morgan.chrom)){
          est.conversion.function <- 1/bp.per.Morgan.chrom[chrom.i]
        }
      }
      
      chromosome.conversion.factor <- c(chromosome.conversion.factor, 1/bp.per.Morgan)
      
      # bp.per.Morgan
      # 
      # bp.per.Morgan <- 1000
      # est.conversion.function <- function(bp, bpm = bp.per.Morgan){bp / bpm}
      # est.conversion.function(10)
      # 
      # bp <- 10
      # 
      # 
      # predict()
      
      #cat(paste0("Skipped chromosome ", chrom.i, " because not enough data available.\n"))
      cat(paste0("Used bp conversion factor for chromosome ", chrom.i, " because not enough data available.\n"))
      
    } else{
      
      # https://stackoverflow.com/questions/25447999/how-to-make-monotonic-increasing-smooth-spline-with-smooth-spline-function
      # ? stats::splinefun
      
      bp <- as.numeric(map.chrom[,4])
      map.chrom <- map.chrom[order(bp),]
      M <- as.numeric(map.chrom[,3])
      bp <- as.numeric(map.chrom[,4])
      
      est.conversion.function <- stats::splinefun(x = bp, y = M
                                                  #, method = "monoH.FC"
                                                  , method = "hyman"
                                                  , ties = list("ordered", mean)
      )
      
      plot(bp, M, xlab = "physical position in bp", ylab = "genetic position in Morgan",
           main = paste0("Chromosome ", chrom.i))
      lines(bp, est.conversion.function(bp), col = "red")
      legend("topleft", legend = c("estimated spline function"), 
             lwd = c(1), 
             col = "red")
      
      if(estimate.conversion.factor & length(bp.per.Morgan.chrom) <= 1){
        chromosome.conversion.factor <- c( chromosome.conversion.factor, 1 / ( (max(bp) - min(bp)) / (max(M) - min(M)) ) )
      } else{
        chromosome.conversion.factor <- c(chromosome.conversion.factor, 1/bp.per.Morgan)
      }
      
    }
    
    cleaned.map <- rbind(cleaned.map, map.chrom)
    list.conversion.functions.chromosomes[[chrom.i]] <- est.conversion.function
  }
  
  names(list.conversion.functions.chromosomes) <- unique(map[,1])
  
  if(!is.null(n.chromosomes) & n.chromosomes > length(unique(map[,1])) ){
    # here, we just insert the bp to cM conversion factor for every chromosome that is not on the map
    est.conversion.function <- 1/bp.per.Morgan
    for(chrom.i in (length(unique(map[,1])) +1):n.chromosomes ){
      list.conversion.functions.chromosomes[[chrom.i]] <- est.conversion.function
      chromosome.conversion.factor <- c(chromosome.conversion.factor, 1/bp.per.Morgan)
    }
    
    names(list.conversion.functions.chromosomes) <- c(
      names(list.conversion.functions.chromosomes)[1:length(unique(map[,1]))],
      (length(unique(map[,1])) +1):n.chromosomes
    )
    
  }
  
  if(return.cleaned.map){
    return(cleaned.map)
  }
  return(
    list(list.conversion.functions.chromosomes = list.conversion.functions.chromosomes,
         list.min.max.bp.chromosomes = list.min.max.bp.chromosomes,
         general.conversion.factor = 1/bp.per.Morgan,
         chromosome.conversion.factor = chromosome.conversion.factor)
  )
}


# this function takes physical positions and assigns genetic (Morgan) positions
# to them based on an estimated function (i.e., a map with SNP with physical and genetic coordinates)
# Also, it can deal with the function just being a conversion factor
assign.genetic.positions <- function(output.function.bp.Morgan, map){
  map <- map[order(as.numeric(map[,4])),]
  map <- map[order(as.numeric(map[,1])),]
  
  unq.chrom.map <- unique(map[,1])
  unq.chrom.output <- unique(names(output.function.bp.Morgan[[1]]))
  
  map.out <- c()
  
  chrom.i <- 2
  for(chrom.i in 1:length(unq.chrom.map)){
    
    if(unq.chrom.map[chrom.i] %in% unq.chrom.output){
      bp <- as.integer( map[map[,1] == unq.chrom.map[chrom.i],4] )
      
      IN <- output.function.bp.Morgan[[1]][[which(unq.chrom.output == unq.chrom.map[chrom.i])]]
      
      # if there was not enough data, then the output of function.bp.Morgan is just a number
      # which is the conversion factor
      is.conversion.factor <- length(IN) == 1 & length(bp) > 1 & class(IN) != "function"
      # class(IN)
      
      if(is.conversion.factor){
        Morgan <- bp * IN
      } else{
        
        min.max.template.gen.map <- output.function.bp.Morgan[[2]][[which(unq.chrom.output == unq.chrom.map[chrom.i])]]
        on.template.gen.map <- bp >= min(min.max.template.gen.map) & bp <= max(min.max.template.gen.map)
        
        any.before <- any(bp < min(min.max.template.gen.map))
        any.after <- any(bp > max(min.max.template.gen.map))
        
        any.before
        any.after
        
        if(any.before == FALSE & any.after == FALSE){
          Morgan <- IN(x = bp)
          
        } else {
          Morgan.on.template.gen.map <- IN(x = bp[on.template.gen.map])
          Morgan.on.template.gen.map.org <- Morgan.on.template.gen.map
          
          if(any.before){
            length(Morgan.on.template.gen.map)
            Morgan.on.template.gen.map <- Morgan.on.template.gen.map - min(Morgan.on.template.gen.map) 
          }
          
          Morgan.before.bp.on.template.gen.map <- c()
          if(any.before){
            Morgan.before.bp.on.template.gen.map <- 
               (bp[bp < min(min.max.template.gen.map)] - min(bp[bp < min(min.max.template.gen.map)]) ) * output.function.bp.Morgan[[4]][[chrom.i]]
            
            x.before.to.on.template <- c(
              min(bp[on.template.gen.map]), 
              min(bp[on.template.gen.map]) + ( min(bp[on.template.gen.map]) - max(bp[bp < min(min.max.template.gen.map)]) )
            )
            M.between.before.and.on.template <- 0
            M.between.before.and.on.template <- diff(IN(x = x.before.to.on.template)) 
          }
          
          Morgan.after.bp.on.template.gen.map <- c()
          if(any.after){
            Morgan.after.bp.on.template.gen.map <- 
              (bp[bp > max(min.max.template.gen.map)] - min(bp[bp > max(min.max.template.gen.map)]) ) * output.function.bp.Morgan[[4]][[chrom.i]]
            
            x.after.to.on.template <- c(
              max(bp[on.template.gen.map]), 
              max(bp[on.template.gen.map]) + ( min(bp[bp > max(min.max.template.gen.map)]) - max(bp[on.template.gen.map]) )
            )
            M.between.after.and.on.template <- 0
            M.between.after.and.on.template <- diff(IN(x = x.after.to.on.template)) 
          }
          
          
          
          if(any.before){
            Morgan.on.template.gen.map <- 
              max(Morgan.before.bp.on.template.gen.map) + M.between.before.and.on.template + Morgan.on.template.gen.map
          }
          
          if(any.after){
            Morgan.after.bp.on.template.gen.map <- 
              max(Morgan.on.template.gen.map) + M.between.after.and.on.template + Morgan.after.bp.on.template.gen.map
          }
          
          Morgan <- c(Morgan.before.bp.on.template.gen.map, Morgan.on.template.gen.map, Morgan.after.bp.on.template.gen.map)
          
          length(on.template.gen.map)
          length(Morgan)
          
          
          first.snp.on.template <- min(which(on.template.gen.map))
          
          # so if the first SNP moved forward, then we readjust the position
          # it may also be that there are so many SNPs in front of the first SNP on teh template,
          # that the template SNP moved behind. In that case, we cannot adjust because it would mean 
          # that the SNPs in front would get a negative position
          if(Morgan[first.snp.on.template] < Morgan.on.template.gen.map.org[1]){
            Morgan <- Morgan + diff(c(Morgan[first.snp.on.template], Morgan.on.template.gen.map.org[1]))
          } else{
            cat(paste0("Genetic positions of chromosome ", chrom.i, " are shifted by ", 
                       round(diff(c(Morgan[first.snp.on.template], Morgan.on.template.gen.map.org[1]))*-1,6),
                       " Morgan.\n"))
          }
        }
      }
      
      # if for some reason the first position is negative, then we need to adjust it so that it 
      # becomes 0
      if(min(Morgan) < 0){
        Morgan <- Morgan + min(Morgan)*-1
      }
      
      # I think MoBPS does not accept genetic positions to be 0
      # so that is why we need to add a small step on top
      if(min(Morgan) == 0){
        d <- diff(Morgan)
        Morgan <- Morgan + min(d[d != 0])
      }
      
      out.map.chrom <- cbind(
        Chromosome = map[map[,1] == unq.chrom.map[chrom.i],1], 
        `SNP ID` = map[map[,1] == unq.chrom.map[chrom.i],2],
        `Morgan position` = Morgan,
        `BP position` = bp
      )
      
    } else{
      next
    }
    
    map.out <- rbind(map.out, out.map.chrom)
    
    cat(paste0("Finished chromosome ", chrom.i, ".\n"))
  }
  
  return(map.out)
}


# function.bp.Morgan() also does some primitive clean up of the genetic map
# that is used as a template. For example, it removes outliers and SNPs that have rank issues 
map.groenen.cleaned <- function.bp.Morgan(map = map.groenen, bp.per.Morgan = 30000000, n.chromosomes = 30, return.cleaned.map=TRUE)

estimated.coversion.function <- function.bp.Morgan(map = map.groenen, bp.per.Morgan = 30000000, 
                                                   n.chromosomes = 30, return.cleaned.map=FALSE, 
                                                   estimate.conversion.factor = TRUE)

# MoBPSmaps::map_chicken1 is the Affymetrix Chicken600K Array
# we only have physical positions but no genetic positions (Morgan)
map.IN <- assign.genetic.positions(output.function.bp.Morgan = estimated.coversion.function, map = MoBPSmaps::map_chicken1)

colros <- rep("grey", 28)
colros[estimated.coversion.function[[4]] == estimated.coversion.function$general.conversion.factor] <- "green"

barplot(height = (1/estimated.coversion.function[[4]])[1:28], 
        names.arg = 1:28, 
        xlab = "chromosome", 
        ylab = "bp per Morgan"
        , ylim = c(0, 50000000 )#max(1/e[[4]])) )
        , col = colros
)
#abline(h = 1/e$general.conversion.factor, col = "green", lwd = 2, lty = "dashed")
legend("topright"
       #, col = c("grey", "green")
       , legend = c("chromosome-specific conversion factor", "genome-wide conversion factor")
       , fill = c("grey", "green"))


map.IN2 <- assign.genetic.positions(output.function.bp.Morgan = estimated.coversion.function, 
                                   map = rbind(map.groenen, MoBPSmaps::map_chicken1[,1:4])
                                     )

# we do this so that the coordinates are the same (they might have gotten shifted during assignment)
map.groenen.cleaned2 <- map.IN2[map.IN2[,2] %in% map.groenen.cleaned[,2],]

map.groenen.cleaned <- map.groenen.cleaned2
map.IN <- map.IN2[!(map.IN2[,2] %in% map.groenen.cleaned2[,2]),]

i <- 1
for(i in 1:28){
  chromosome <- i
  if(is.numeric(estimated.coversion.function[[1]][[chromosome]]) == 1){next}
  
  pos.groenen <- map.groenen.cleaned[,1] == chromosome
  bp <- map.groenen.cleaned[pos.groenen,4]
  M <- map.groenen.cleaned[pos.groenen,3]
  
  pos <- map.IN[,1] == chromosome
  
  plot(x = bp, 
       y = M, 
       main = paste0("Chromosome ", chromosome),
       ylab = "genetic position in Morgan",
       xlab = "physical position in bp"
       , xlim = c(min(c(as.numeric(map.IN[pos,4]), as.numeric(map.groenen.cleaned[pos.groenen,4]))), 
                  max(c(as.numeric(map.IN[pos,4]), as.numeric(map.groenen.cleaned[pos.groenen,4]))))
       , ylim = c(min(c(as.numeric(map.IN[pos,3]), as.numeric(map.groenen.cleaned[pos.groenen,3]))), 
                  max(c(as.numeric(map.IN[pos,3]), as.numeric(map.groenen.cleaned[pos.groenen,3]))))
  )
  
  lines(bp, estimated.coversion.function[[1]][[chromosome]](bp), col = "red")
  legend("topleft", legend = c("estimated spline function"), 
         lwd = c(1), 
         col = "red"
         , bg = "beige")
  
  pos <- map.IN[,1] == chromosome
  map.chrom <- map.IN[pos,]
  bp <- as.integer(map.chrom[,4])
  Morgan <- as.numeric(map.chrom[,3])
  
  r <- ceiling(seq(from = 1, to = length(bp), length.out = min(c(150, length(bp)))))
  points(x = bp[r], y = Morgan[r], pch = 4, col = "dodgerblue", lwd = 2)
  legend("bottomright", pch = 4, col = "dodgerblue", legend = 'Affymetrix Chicken600K Array \ngenetic positions estimated\n', lwd = 2, lty = 0
         , bg = "beige")
  
  Sys.sleep(0.2)
}

# sometimes, you can see deviates from the line lien to the blue points. That is no problem. this is because some
# distance was added to avoid having some positions as negative coordinates
#
# Manual curation is still adviced. Positions at chromsome ends are estimated based on whole chromosome conversion factor
# also, some SNPs seem to be off (see chromosome 22)