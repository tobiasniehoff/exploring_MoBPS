# functions for OGC script
# there is some overlap with function for them helper function file


# this is a primitive implementation for OCG based on optiSel
# examples for code can be found in Robin Wellmann's publication 2019 zu optiSel
#
# ub.fPED  <- av.kin  + (1-av.kin)/(2*Ne*generation.interval)
opt.contr <- function(population, database, kinship = NULL, pop.size = NULL, bvs,
                      #ub.kin=NULL, lb.BV=NULL, ub.BV=NULL, 
                      objective = "max.EBV", verbose = TRUE,
                      con=NULL, solver = "cccp"# use"slsqp" if cccp doesn't work
){
  library("optiSel")
  library("data.table")
  if(is.null(con)){
    stop("Need to provide list with contraints.\n
         Look in optiSel publication.\n")
  }
  
  if(is.null(kinship)){
    kinship <- kinship.exp(population, database = database)
  }
  kinship <- qc.combination.matrix(kinship, cut.symmetric.part = F)
  ids <- get.id(population, database = database, use.id = T)
  
  sex <- get.sex(population, database = database)
  sex[sex == 1] <- "male"
  sex[sex == 2] <- "female"
  
  if(is.null(pop.size)){
    pop.size <- 1000
  }
  
  phen <- data.table(Indiv = ids, Breed = "BreedA", 
                     Sex = sex, EBV = as.vector(bvs), 
                     isCandidate = TRUE)
  cand <- candes(phen = phen, fPED = kinship, 
                 N = pop.size, 
                 quiet = !verbose)
  
  fit <- opticont(objective, cand = cand, con = con, solver=solver, quiet = !verbose)
  
  return(fit)
}

# db <- get.database(pop, gen = length(pop$breeding))
# bvs <- get.bv(pop, database = db, use.id = T)[1,]
# kin.exp <- kinship.exp(pop, database = db)
# mean(kin.exp[upper.tri(kin.exp, diag = F)])
# 
# # gives different results depending on the population size
# # cand <- candes(phen = phe, fPED = kin.exp, N = pop.size)
# # cand$mean$fPED
# 
# max.kin.next.gen <- mean(kin.exp[upper.tri(kin.exp, diag = F)]) + 
#   (1-mean(kin.exp[upper.tri(kin.exp, diag = F)]))*0.01
# max.kin.next.gen <- mean(kin.exp[upper.tri(kin.exp, diag = F)]) + 0.05
# 
# con <- list(ub.fPED = max.kin.next.gen)
# 
# res <- opt.contr(pop, database = db, kinship = kin.exp, 
#                  objective = "max.EBV", con = con, 
#                  bvs = bvs, solver = "cccp")
# 
# 
# res$info
# 
# res$mean
# Candidate <- res$parent
# Candidate[Candidate$oc>0.01,c("Sex","EBV","oc")]

# quality control of combination matrix. Expects integers as IDs
qc.combination.matrix <- function(combination.matrix, reorder = TRUE, cut.symmetric.part = TRUE){
  parents.id <- sort(as.integer(unique(c(rownames(combination.matrix), colnames(combination.matrix)))))
  
  if(!any(rownames(combination.matrix) %in% colnames(combination.matrix))){
    ncmbm <- matrix(0, nrow = length(parents.id), ncol = length(parents.id))
    ncmbm[1:nrow(combination.matrix), 1:ncol(combination.matrix)] <- combination.matrix
    ncmbm[(ncol(combination.matrix)+1):ncol(ncmbm), 
          (nrow(combination.matrix)+1):nrow(ncmbm)] <- t(combination.matrix)
    
    rownames(ncmbm) <- c(rownames(combination.matrix), colnames(combination.matrix))
    colnames(ncmbm) <- c(colnames(combination.matrix), rownames(combination.matrix))
    combination.matrix <- ncmbm
  }
  
  if(reorder){
    # reordering of matrix
    combination.matrix <- combination.matrix[,match(colnames(combination.matrix),parents.id)]
    combination.matrix <- combination.matrix[match(rownames(combination.matrix),parents.id),]
  }
  
  # if matrix is supposed to be symeteric (DMaxSire is same as SirexDam), 
  # then you can cut out half of teh matrix
  if(cut.symmetric.part){
    combination.matrix[(combination.matrix + t(combination.matrix)) != 0] <- 1
    combination.matrix[lower.tri(combination.matrix)] <- 0
  }
  
  return(combination.matrix)
}

# THIS IS THE PRIMIIVE ALGORITHM FOR CONSIDERING NUMBER OF PARENTS AND HAVING NICE SPEED
#
# to be sure only the right individuals are removed, you need to remove
# one animal, recalculate oc, remove one again and so on.
# But that takes very long and many animals have only a tiny contribution
# This function removes all animals that have tiny contributions.
# For that, it sorts the aniamls per sex from most contributing to least.
# Then, the cumulative sum is calculated and only those animals are kept
# that have a cumulative sum lower the threshold
# The threshold is per sex. A male can at most only contribute 50 percent
# but here this is sen as 100% (oc *2)
# The motivation for this function is that many animals only have atiny contribution
# use a threshold value of 0.9999 or so.
# n.min.animal 2 means at least twice as many animals of the sex must be left over
# to use this strategy. Otherwise, nothing is removed
# max.rm.animal means how many animals should be removed at once at maximum
remove.perc.ind.lowest.contr <- function(result.opt.contr, n.sel.males, 
                                         n.sel.females, ids, 
                                         thres.min.cumsum.contr = 1,
                                         n.min.animal = c(2.5,1.5),
                                         max.rm.animal = 10){
  res <- result.opt.contr
  res$parent <- res$parent[res$parent[,"Indiv"] %in% ids,]
  
  # this is checking whether it is possible to remove individuals
  if(((length(res$parent[res$parent$Sex == "male", "oc"]) > n.sel.males) + 
      (length(res$parent[res$parent$Sex == "female", "oc"]) > n.sel.females)) > 0 |
     thres.min.cumsum.contr < 1){
    
    pos <- order(res$parent[res$parent$Sex == "male", c("oc")], 
                 decreasing = T)
    
    repos.m <- res$parent[res$parent$Sex == "male", c("Indiv", "oc")][pos,]
    CUMSUM.m <- cumsum(repos.m[,"oc"])*2
    CUMSUM.m[CUMSUM.m > 1] <- 1
    
    pos <- order(res$parent[res$parent$Sex == "female", c("oc")], 
                 decreasing = T)
    repos.f <- res$parent[res$parent$Sex == "female", c("Indiv", "oc")][pos,]
    CUMSUM.f <- cumsum(repos.f[,"oc"])*2
    CUMSUM.f[CUMSUM.f > 1] <- 1
    
    if(sum(c(CUMSUM.m, CUMSUM.f) > thres.min.cumsum.contr) > max.rm.animal){
      
      names(CUMSUM.m) <- repos.m[, "Indiv"]
      names(CUMSUM.f) <- repos.f[, "Indiv"]
      
      CUMSUM <- c(CUMSUM.m, CUMSUM.f)
      rm.ids <- names(
        tail(which(CUMSUM[order(CUMSUM, decreasing = F)] >= thres.min.cumsum.contr), 
             n = max.rm.animal))
      
      ids.m <- names(CUMSUM.m)[!names(CUMSUM.m) %in% rm.ids]
      ids.f <- names(CUMSUM.f)[!names(CUMSUM.f) %in% rm.ids]
    } else{
      ids.m <- repos.m[CUMSUM.m <= thres.min.cumsum.contr,][,"Indiv"]
      ids.f <- repos.f[CUMSUM.f <= thres.min.cumsum.contr,][,"Indiv"]
    }
    
    # this is checking the n.min.animals condition
    if(length(ids.m) <= (n.sel.males * n.min.animal[1]) |
       length(ids.f) <= (n.sel.females * n.min.animal[2])){
      out.ids <- ids
    } else{
      out.ids <- ids[ids %in% c(ids.m, ids.f)]
    }
    
    #length(ids.f) + length(ids.m)
    # here I reorder the selected animals so that the order matches the input order
    
    #length(out.ids)
    
    cat(paste0("Removed ", length(ids)-length(out.ids), " animals with shortcut.\n"))
    ids <- out.ids
  }
  
  return(ids)
}
# a <- remove.perc.ind.lowest.contr(res, n.sel.males = 20, 
#                                   n.sel.females = 30, ids = ids, 
#                                   thres.min.cumsum.contr = 0.99999)
# 
# length(a)  
# a
# w <- res$parent[res$parent$Sex == "female",c("Indiv", "oc")]
# 
# q <- w[as.integer(w[,"Indiv"]) %in% a,]
# 
# sum(q[,"oc"]*2)

# this function removes the individual with the lowest contribution
# given the minimum number of to be selected males and females is not already reached
remove.ind.lowest.contr <- function(result.opt.contr, n.sel.males, n.sel.females, ids){
  
  res <- result.opt.contr
  res$parent <- res$parent[res$parent[,"Indiv"] %in% ids,]
  male.smaller <- min(res$parent[res$parent$Sex == "male", "oc"]) < min(res$parent[res$parent$Sex == "female", "oc"])
  
  if(((length(res$parent[res$parent$Sex == "male", "oc"]) > n.sel.males) + 
      (length(res$parent[res$parent$Sex == "female", "oc"]) > n.sel.females)) > 0){
    
    if(male.smaller){
      cat("A male has the lowest contribution.\n")
      # if the list is not already reduced to the final number of males
      if(length(res$parent[res$parent$Sex == "male", "oc"]) > n.sel.males){
        # remove lowest male
        rm.ind <- res$parent[res$parent$Sex == "male",][
          which.min(res$parent[res$parent$Sex == "male", "oc"]), "Indiv"]
        cat("Removed a male.\n")
      } else{
        # this else is for the case that a male has lowest contribution but minimum 
        # number of selected males is already reached
        # so lowest female needs to be removed
        rm.ind <- res$parent[res$parent$Sex == "female",][
          which.min(res$parent[res$parent$Sex == "female", "oc"]), "Indiv"]
        cat("Removed a female.\n")
      }
    } else{
      cat("A female has the lowest contribution.\n")
      # this is for the case that a female contribution is smaller
      # if the list is not already reduced to the final number of females
      if(length(res$parent[res$parent$Sex == "female", "oc"]) > n.sel.females){
        # remove lowest female
        rm.ind <- res$parent[res$parent$Sex == "female",][
          which.min(res$parent[res$parent$Sex == "female", "oc"]), "Indiv"]
        cat("Removed a female.\n")
      } else{
        # this else is for the case that a female has lowest contribution but minimum 
        # number of selected females is already reached
        # so lowest male needs to be removed
        rm.ind <- res$parent[res$parent$Sex == "male",][
          which.min(res$parent[res$parent$Sex == "male", "oc"]), "Indiv"]
        cat("Removed a male.\n")
      }
    }
    cat(paste0("Removed individual ", rm.ind, "\n"))
    ids <- ids[-which(ids == rm.ind)]
  }
  
  return(ids)
}
# n.sel.males <- 5
# n.sel.females <- 20
# res <- opt.contr(pop, database = db, kinship = kin.exp, objective = "max.EBV", 
#                  solver = "cccp", bvs = bvs, con = con)
# ids <- get.id(pop, database = db, use.id = T)

# res is the object given by opticont
get.vector.occurence <- function(res, population.size){
  Sire.vec <- 
    rep(sel.m[,"Indiv"], times = round((sel.m[,"oc"]/sum(sel.m[,"oc"]))*population.size))
  resample <- function(x, ...) {x[sample.int(length(x), ...)]} # idea from ?sample
  Sire.vec <- resample(Sire.vec)
  if(length(Sire.vec) < population.size){
    
    # fill in the missing positions randomly with probabilities from oc
    Sire.vec <- c(Sire.vec,
                  sample(as.integer(sel.m[,"Indiv"]), 
                         prob = sel.m[,"oc"]/sum(sel.m[,"oc"]), 
                         replace = T, 
                         size = population.size-length(Sire.vec)))
  }
  if(length(Sire.vec) > population.size){
    Sire.vec <- Sire.vec[-sample(1:population.size, size = length(Sire.vec) - population.size)]
  }
  return(Sire.vec)
}
