

library(MoBPS) # 1.9.23
library(stringr)
library(readxl)
library("optiSel")
library("data.table")
library(writexl)


odir <- getwd()
#(odir <- dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(odir)
getwd()

source("functions_for_OGC_script.R")
setwd(odir)
################################################################################
args_script <- commandArgs(trailingOnly = TRUE)
index_for_scenario <- as.numeric(args_script[1])

rowi <- index_for_scenario
rowi <- as.numeric(Sys.getenv(c("SLURM_ARRAY_TASK_ID")))
ncores <- as.numeric(Sys.getenv(c("SLURM_CPUS_PER_TASK")))

name_parameter_table <- args_script[2]
# name_parameter_table <- 'OGC_styles_simple_investigation1'
# rowi <- 1
parameters <- as.data.frame(read_excel(paste0(name_parameter_table,".xlsx"), skip = 1))

################################################################################
if(!dir.exists(name_parameter_table)){
  dir.create(name_parameter_table)
}
setwd(paste0('./',name_parameter_table))
getwd()

n.sel.males <- parameters$n.sel.males[rowi]
n.sel.females <- parameters$n.sel.females[rowi]
population.size <- parameters$population.size[rowi]
set.seed(parameters$seed[rowi])
OGC.type <- parameters$OGC.type[rowi]
generations <- parameters$generations[rowi]
basis <- parameters$basis[rowi]
projection <- parameters$projection[rowi]
ogc.relationship.matrix <- parameters$ogc.relationship.matrix[rowi]
one.male.per.female <- parameters$one.male.per.female[rowi]
thresh.frac.shortcut <- as.numeric(parameters$thresh.frac.shortcut[rowi])
gen.gain <- parameters$gen.gain[rowi]
objective <- parameters$objective[rowi]
use.invalid.solution <- parameters$use.invalid.solution[rowi]
Ne <- parameters$Ne[rowi]
time.sim.start <- Sys.time()

pop <- creating.diploid(chr.nr = 30, 
                        n.additive = 5000,
                        nindi = 1000, 
                        nsnp = 55000,
                        chromosome.length = 1,
                        snps.equidistant = T,
                        base.bv = 0,
                        var.target = 10
)

# generate some relationship
for(i in 1:10){
  pop <- breeding.diploid(pop, 
                          selection.size = c(100,100), 
                          breeding.size = c(500,500), 
                          delete.same.origin = T, verbose = F)
}


# drive population
for(i in 1:5){
  pop <- breeding.diploid(pop, 
                          selection.size = c(n.sel.males,n.sel.females), 
                          breeding.size = c(population.size/2,population.size/2), 
                          delete.same.origin = T, 
                          selection.criteria = c("bv", "bv"), verbose = F
                          #,generation.cores = 3
                          )
}

#2^15

pop_save <- pop
pop <- pop_save

get.effective.size(pop, gen = length(pop$breeding))

set.seed(parameters$seed[rowi])
L <- 1

(k <- kinship.emp.fast(pop, gen = length(pop$breeding), ibd.obs = 1000, hbd.obs = population.size))
vec.av.kin <- k[1]
vec.av.inbreeding <- k[2]
vec.av.bv <- mean(get.bv(pop, gen = length(pop$breeding))[1,])

pop2 <- breeding.diploid(pop, 
                         breeding.size = 500, 
                         selection.m.database = get.database(pop, gen = length(pop$breeding))[1,],
                         selection.f.database = get.database(pop, gen = length(pop$breeding))[2,], 
                         max.offspring = ceiling(c(500/(population.size/2), 500/(population.size/2))), 
                         max.mating.pair = 1, verbose = F)

vec.gen.var <- var(get.bv(pop2, gen = length(pop2$breeding))[1,])

valid.solution <- c()
#i <- 1

for(i in 1:generations){
  
  calculating <- FALSE
  if(projection == "every.generation"){
    calculating <- TRUE
  }
  if(projection == "linear.projection" & i == 1){
    calculating <- TRUE
  }
  if(calculating){
    if(basis == "pop.av"){
      # the below calculation is more accurate
      sex <- get.sex(pop, gen = length(pop$breeding))
      sex[sex == 1] <- "male"
      sex[sex == 2] <- "female"
      
      ids <- as.integer(get.id(pop, gen = length(pop$breeding), use.id = T))
      
      # creating rel matrix
      # bve rel matrix is not identical to kinship matrix
      # obtained with kinship.exp kinship.exp
      if(ogc.relationship.matrix == "pedigree"){
        # if we need a pedigree rel matrix, this is faster
        kin.exp <- kinship.exp(pop, gen = length(pop$breeding))
      } else{
        kin.exp <- breeding.diploid(pop, 
                                    export.relationship.matrix = T, 
                                    bve = T, 
                                    bve.gen = length(pop$breeding), 
                                    relationship.matrix = ogc.relationship.matrix,
                                    remove.effect.position = T)
      }
      row.names(kin.exp) <- sort(ids)
      colnames(kin.exp) <- sort(ids)
      
      # mean(kin.exp[upper.tri(kin.exp)])
      # ke <- kinship.exp(pop, gen = length(pop$breeding))
      # mean(ke[upper.tri(ke)])
      
      # max.kin.next.gen <- mean(kin.exp[upper.tri(kin.exp, diag = T)]) + 0.01
      bvs <- get.bv(pop, gen = length(pop$breeding), use.id = T)[1,]
      gen.sd <- sd(bvs)
      phen <- data.table(Indiv = ids, Sex = sex, EBV = as.vector(bvs), isCandidate = TRUE)
      cand <- candes(phen = phen, fPED = kin.exp, quiet = T)
      
      # this is the average kinship in the population in the next generation if 
      # animals were mated at random with males 50% contribution and females 50% contribution
      # specifying the population size accounts for drift
      cand$mean$fPED
    }
    if(basis == "sel.animals.av"){
      sex <- get.sex(pop, gen = length(pop$breeding), use.id = T)
      bvs <- get.bv(pop, gen = length(pop$breeding), use.id = T)
      
      ids.sel.m <- sort(as.integer(names(bvs[,sex == 1])[order(bvs[,sex == 1], decreasing = T)[1:n.sel.males]]))
      db.sel.m <- get.database(pop, id = ids.sel.m)
      
      ids.sel.f <- sort(as.integer(names(bvs[,sex == 2])[order(bvs[,sex == 2], decreasing = T)[1:n.sel.females]]))
      db.sel.f <- get.database(pop, id = ids.sel.f)
      
      db <- rbind(db.sel.m, db.sel.f)
      ids <- as.integer(get.id(pop, database = db, use.id = T))
      
      if(ogc.relationship.matrix == "pedigree"){
        # if we need a pedigree rel matrix, this is faster
        kin.exp <- kinship.exp(pop, database = db)
      } else{
        kin.exp <- breeding.diploid(pop, 
                                    export.relationship.matrix = T, 
                                    bve = T, 
                                    bve.database = db, 
                                    relationship.matrix = ogc.relationship.matrix,
                                    remove.effect.position = T)
      }
      row.names(kin.exp) <- sort(ids)
      colnames(kin.exp) <- sort(ids)
      
      sex <- get.sex(pop, database = db, use.id = T)
      sex[sex == 1] <- "male"
      sex[sex == 2] <- "female"
      
      bvs <- get.bv(pop, database = db, use.id = T)[1,]
      gen.sd <- sd(bvs)
      phen <- data.table(Indiv = ids, Sex = sex, EBV = as.vector(bvs), isCandidate = TRUE)
      
      cand <- candes(phen = phen, fPED = kin.exp, quiet = F)
      #get.effective.size(pop, database = db)
      cand$mean$fPED
    }
  }
  
  if(projection == "linear.projection"){
    #max.kin.next.gen <- cand$mean$fPED + 0.01*i
    max.kin.next.gen <- 1-(1-cand$mean$fPED) * (1- (1/(2*Ne))) ^ (((1+i)-1)/L)
    
    min.bv.next.gen <- cand$mean$EBV + gen.sd*gen.gain*i
  }else{
    #max.kin.next.gen <- cand$mean$fPED + 0.01
    max.kin.next.gen <- 1-(1-cand$mean$fPED) * (1- (1/(2*Ne))) ^ ((2-1)/L)
    # (max.kin.next.gen.n <- pop.av.kin + (1-pop.av.kin)/(2*Ne))
    # this gives same result if L is 1. Easier to understand
    min.bv.next.gen <- cand$mean$EBV + gen.sd*gen.gain
  }
  
  if(objective == "max.EBV"){
    con <- list(ub.fPED = max.kin.next.gen)
  } else{
    con <- list(lb.EBV = min.bv.next.gen)
  }
  
  if(one.male.per.female){
    offspring.per.mating <- population.size/n.sel.females
    # defining upper bounds for contribution for females
    ub <- setNames(
      rep(offspring.per.mating/(2*population.size), 
          population.size/2),
      get.id(pop, 
             gen = length(pop$breeding), 
             use.id = T)[get.sex(pop, 
                                 gen = length(pop$breeding)) == 2]
    )
    if(objective == "max.EBV"){
      con <- list(ub.fPED = max.kin.next.gen, ub = ub)
    } else{
      con <- list(lb.EBV = min.bv.next.gen, ub = ub)
    }
  }
  
  db <- get.database(pop, gen = length(pop$breeding))
  ids <- get.id(pop, database = db, use.id = T)
  
  if(OGC.type == "OGC.1type"){
    # OGC type 1 is selecting the best animals and then restricting inbreeding
    sex <- get.sex(pop, database = db)
    bvs <- get.bv(pop, database = db, use.id = T)
    
    ids.sel.m <- sort(as.integer(names(bvs[,sex == 1])[order(bvs[,sex == 1], decreasing = T)[1:n.sel.males]]))
    db.sel.m <- get.database(pop, id = ids.sel.m)
    
    ids.sel.f <- sort(as.integer(names(bvs[,sex == 2])[order(bvs[,sex == 2], decreasing = T)[1:n.sel.females]]))
    db.sel.f <- get.database(pop, id = ids.sel.f)
    
    db <- rbind(db.sel.m, db.sel.f)
    ids <- as.integer(get.id(pop, database = db, use.id = T))
    if(ogc.relationship.matrix == "pedigree"){
      # if we need a pedigree rel matrix, this is faster
      kin.exp <- kinship.exp(pop, database = db)
    } else{
      kin.exp <- breeding.diploid(pop, 
                                  export.relationship.matrix = T, 
                                  bve = T, 
                                  bve.database = db, 
                                  relationship.matrix = ogc.relationship.matrix,
                                  remove.effect.position = T)
    }
    row.names(kin.exp) <- sort(ids)
    colnames(kin.exp) <- sort(ids)
    
    bvs <- get.bv(pop, database = db, use.id = T)[1,]
    #get.effective.size(pop, database = db)
    # sex <- get.sex(pop, database = db)
    # sex[sex == 1] <- "male"
    # sex[sex == 2] <- "female"
    # 
    # phen <- data.table(Indiv = get.id(pop, database = db, use.id = T),
    #                    Sex = sex, isCandidate = TRUE)
    # cand <- candes(phen = phen, fPED = kin.exp, quiet = F)
    if(one.male.per.female){con$ub <- con$ub[names(con$ub) %in% ids.sel.f]}
    # NOTE: Sometimes, the result is not valid 
    # this is likely because the effective population
    # size is not large enough
    #
    # you can check the effective size with
    #get.effective.size(pop, database = db)
    #
    # Anyway, if there is no valid solution, there is no guarantee
    # the the output solution is as close to the optimal solution as possible
    # It could also be complete rubbish
    
    res <- opt.contr(pop, database = db, kinship = kin.exp, objective = objective, 
                     solver = "slsqp", bvs = bvs, con = con,
                     verbose = T)
    
    # cm <- sum(res$parent[res$parent$Sex == "male","oc"])
    # cf <- sum(res$parent[res$parent$Sex == "female","oc"])
    # cm/(cm+cf)
    # cf/(cm+cf)
    
    # if slsqp fails, its solution is worse than when cccp2 fails
    if(!res$info$valid){
      cat("Solver slsqp failed. Using cccp2 now. \n")
      res.cccp2 <- try(opt.contr(pop, database = db, kinship = kin.exp, objective = objective, 
                         solver = "cccp2", bvs = bvs, con = con, pop.size = population.size,
                         verbose = T), silent = T)
      # so, if cccp2 is better, use solution from cccp2
      if(res.cccp2$summary[res.cccp2$summary$Var == "fPED", "Val"] < 
         res$summary[res$summary$Var == "fPED", "Val"] | res.cccp2$info$valid){
        res <- res.cccp2
      }
    }
    valid.solution <- c(valid.solution, res$info$valid)
    
    res$parent$oc[res$parent$oc < 0] <- 0
    
    sel.m <- res$parent[res$parent$Sex == "male",]
    Sire <- as.integer(get.vector.occurence(sel.m, population.size = population.size))
    
    sel.f <- res$parent[res$parent$Sex == "female",]
    Dam <- as.integer(get.vector.occurence(sel.f, population.size = population.size))
  }
  if(OGC.type == "OGC.2type"){
    
    db <- get.database(pop, gen = length(pop$breeding))
    ids <- as.integer(get.id(pop, database = db, use.id = T))
    
    if(ogc.relationship.matrix == "pedigree"){
      # if we need a pedigree rel matrix, this is faster
      kin.exp <- kinship.exp(pop, gen = length(pop$breeding))
    } else{
      kin.exp <- breeding.diploid(pop, 
                                  export.relationship.matrix = T, 
                                  bve = T, 
                                  bve.gen = length(pop$breeding), 
                                  relationship.matrix = ogc.relationship.matrix,
                                  remove.effect.position = T)
    }
    row.names(kin.exp) <- sort(ids)
    colnames(kin.exp) <- sort(ids)
    
    bvs <- get.bv(pop, database = db, use.id = T)[1,]
    #get.effective.size(pop, database = db)
    res <- opt.contr(pop, database = db, kinship = kin.exp, objective = objective, 
                     solver = "cccp2", bvs = bvs, con = con, pop.size = population.size,
                     verbose = F)
    valid.solution <- c(valid.solution, res$info$valid)
    res$parent$oc[res$parent$oc < 0] <- 0
    
    res$parent[res$parent$Sex == "male","oc"]
    sel.m <- res$parent[res$parent$Sex == "male", ][
      order(res$parent[res$parent$Sex == "male","oc"], 
            decreasing = T)[1:n.sel.males], c("Indiv", "Sex", "oc")]
    Sire <- 
      as.integer(get.vector.occurence(sel.m, population.size = population.size))
    
    sel.f <- res$parent[res$parent$Sex == "female", ][
      order(res$parent[res$parent$Sex == "female","oc"], 
            decreasing = T)[1:n.sel.females], c("Indiv", "Sex", "oc")]
    Dam <- 
      as.integer(get.vector.occurence(sel.f, population.size = population.size))
  }
  if(OGC.type == "OGC.2.2type"){
    # this type is basically two loops
    db <- get.database(pop, gen = length(pop$breeding))
    ids <- as.integer(get.id(pop, database = db, use.id = T))
    all.sol.valid <- c()
    for(l in 1:2){
      db <- get.database(pop, id = ids)
      
      if(ogc.relationship.matrix == "pedigree"){
        # if we need a pedigree rel matrix, this is faster
        kin.exp <- kinship.exp(pop, database = db)
      } else{
        kin.exp <- breeding.diploid(pop, 
                                    export.relationship.matrix = T, 
                                    bve = T, 
                                    bve.database = db, 
                                    relationship.matrix = ogc.relationship.matrix,
                                    remove.effect.position = T)
      }
      row.names(kin.exp) <- sort(ids)
      colnames(kin.exp) <- sort(ids)
      
      bvs <- get.bv(pop, database = db, use.id = T)[1,]
      #get.effective.size(pop, database = db)
      ids.f <- get.id(pop,
                      database = get.database(pop,
                                              gen = length(pop$breeding))[2, ],
                      use.id = T)
      
      ids.sel.f <- ids[ids %in% ids.f]
      if(one.male.per.female){con$ub <- con$ub[names(con$ub) %in% ids.sel.f]}
      res <- opt.contr(pop, database = db, kinship = kin.exp, objective = objective, 
                       solver = "cccp2", bvs = bvs, con = con, pop.size = population.size,
                       verbose = F)
      all.sol.valid <- c(all.sol.valid, res$info$valid)
      
      res$parent$oc[res$parent$oc < 0] <- 0
      
      sel.m <- res$parent[res$parent$Sex == "male", ][
        order(res$parent[res$parent$Sex == "male","oc"], 
              decreasing = T)[1:n.sel.males], c("Indiv", "Sex", "oc")]
      
      sel.f <- res$parent[res$parent$Sex == "female", ][
        order(res$parent[res$parent$Sex == "female","oc"], 
              decreasing = T)[1:n.sel.females], c("Indiv", "Sex", "oc")]
      
      ids <- as.integer(c(sel.m[, "Indiv"], sel.f[, "Indiv"]))
    }
    valid.solution <- c(valid.solution, any(!all.sol.valid))
    Sire <- 
      as.integer(get.vector.occurence(sel.m, population.size = population.size))
    
    Dam <- 
      as.integer(get.vector.occurence(sel.f, population.size = population.size))
    
  }
  if(OGC.type == "OGC.3type"){
    con.store <- con
    con <- con.store
    db <- get.database(pop, gen = length(pop$breeding))
    ids <- as.integer(get.id(pop, database = db, use.id = T))
    all.sol.valid <- c()
    # looping over all animals
    # technically incorrect because some animals are left so it is not ALL
    # but the loop is broken when it doesn't work anymore
    for(j in 1:population.size){
      
      db <- get.database(pop, id = ids)
      
      if(ogc.relationship.matrix == "pedigree"){
        # if we need a pedigree rel matrix, this is faster
        if(j == 1){
          kin.exp <- kinship.exp(pop, database = db)
        }else{
          kin.exp <- kin.exp[,colnames(kin.exp) %in% ids]
          kin.exp <- kin.exp[row.names(kin.exp) %in% ids,]
        }
        
      } else{
        kin.exp <- breeding.diploid(pop, 
                                    export.relationship.matrix = T, 
                                    bve = T, 
                                    bve.database = db, 
                                    relationship.matrix = ogc.relationship.matrix,
                                    remove.effect.position = T)
      }
      row.names(kin.exp) <- sort(ids)
      colnames(kin.exp) <- sort(ids)
      
      bvs <- get.bv(pop, database = db, use.id = T)[1,]
      
      
      ids.sel.f <- ids[ids %in% 
        get.id(pop, gen = length(pop$breeding), use.id = T)[
          get.sex(pop, gen = length(pop$breeding)) == 2]
      ]
      
      if(one.male.per.female){con$ub <- con$ub[names(con$ub) %in% ids.sel.f]}
      # "slsqp"
      t1 <- Sys.time()
      res <- opt.contr(pop, database = db, kinship = kin.exp, objective = objective, 
                       solver = "cccp2", bvs = bvs, con = con, pop.size = population.size, 
                       verbose = F)
      t2 <- Sys.time()
      
      if(j == 1){
        res.store <- data.frame(Indiv = res$parent$Indiv, 
                                Sex = res$parent$Sex,
                                "oc.iter.1" = res$parent$oc)
      }else{
        res.store <- cbind(res.store, 0)
        res.store[res.store[,1] %in% res$parent$Indiv,ncol(res.store)] <- res$parent$oc
        colnames(res.store)[ncol(res.store)] <- paste0("oc.iter",j)
      }
      
      all.sol.valid <- c(all.sol.valid, res$info$valid)
      
      res$parent$oc[res$parent$oc < 0] <- 0
      
      # hist(res$parent$oc)
      # hist(res$parent[res$parent$Sex == "male", c("oc")])
      cumsum(sort(res$parent[res$parent$Sex == "female", c("oc")]*2, decreasing = T)) > thresh.frac.shortcut
      cumsum(sort(res$parent[res$parent$Sex == "male", c("oc")]*2, decreasing = T)) > thresh.frac.shortcut
      
      length(ids)
      ids <- as.integer(res$parent$Indiv)
      new.ids <- ids
      # this removes animals that have a tiny contribution. This is done to speed up the process
      new.ids <- remove.perc.ind.lowest.contr(res, n.sel.males = n.sel.males, 
                                   n.sel.females = n.sel.females, 
                                   ids = ids, 
                                   thres.min.cumsum.contr = thresh.frac.shortcut
                                   )
      length(new.ids)
      
      #get.effective.size(pop, database = get.database(pop, id = res$parent[,"Indiv"]))
      #ld.decay(pop, database = get.database(pop, id = res$parent[,"Indiv"]), plot = T, chromosome = 3)
      
      #res$parent[res$parent$Sex == "male", c("Indiv","oc")][
      #  order(res$parent[res$parent$Sex == "male", "oc"], decreasing = T)[1:5],]
      
      #list.res <- c(list.res, list(res))
      
      # here the individual with the lowest contribution is removed if the minimum 
      # number of to be selected males and females is not yet reached
      # the condition is checking whether new.ids is different from ids.
      # new.ids can only be different, if the function remove.perc.ind.lowest.contr()
      # removed animals
      if(!any(!(ids %in% new.ids))){
        new.ids <- remove.ind.lowest.contr(result.opt.contr = res, 
                                           n.sel.males = n.sel.males, 
                                           n.sel.females = n.sel.females, 
                                           ids = ids)
      }else{
        cat("Used shortcut. \n")
      }
      
      nfinished <- population.size - length(new.ids)
      prog = max(floor(nfinished/(population.size - n.sel.females - n.sel.males)*50),1)
      minutes.to.go <- round(as.numeric(difftime(t2, t1, units = "secs")*
                                          (population.size - n.sel.females - n.sel.males-nfinished))/60,2)
      
      cat('\r', "|", strrep("#",prog),
          strrep(" ", 50-prog), "| ", 2*prog, "%", sep="", 
          " still need ", minutes.to.go, " minutes")
      flush.console()
      cat(" \n\r")
      
      if(identical(new.ids, ids)){
        break
      }
      ids <- new.ids
    }
    #res.shortcut <- res
    #res.all <- res
    #res.store.shortcut <- res.store
    #res.store.all <- res.store
    # res.store
    # res.shortcut$parent$oc
    # res.all$parent$oc
    #identical(res.all, res)
    # res.all$summary
    # res.shortcut$summary
    # res.9$summary
    # res$summary
    #res.9 <- res
    #res.9.store <- res.store
    #res.shortcut$parent$Indiv %in% res.all$parent$Indiv
    
    #length(all.sol.valid)
    valid.solution <- c(valid.solution, !any(!all.sol.valid))
    
    # if invalid solutions should be used and there is an invalid solution,
    # overwrite the current solution with the last valid one
    if(!use.invalid.solution){
      if(any(!all.sol.valid) & any(all.sol.valid)){
        res$parent <- res.store[,c(1,2,max(which(all.sol.valid)))]
        colnames(res$parent)[3] <- 'oc'
      } else{
        res$parent$oc[res$parent$Sex == "male"] <- 2/length(res$parent$Sex == "male")
        res$parent$oc[res$parent$Sex == "female"] <- 2/length(res$parent$Sex == "female")
      }
    }
    
    # # if solution is invalid, all animal get equal contributions
    # if(!use.invalid.solution){
    #   res$parent$oc[res$parent$Sex == "male"] <- 2/length(res$parent$Sex == "male")
    #   res$parent$oc[res$parent$Sex == "female"] <- 2/length(res$parent$Sex == "female")
    # }
    
    # this is comparing the last to the second to last result
    # cor(res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)][1:12],
    # res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)-1][1:12])
    # 
    # res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)][1:12]/
    #   sum(res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)][1:12])
    # 
    # res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)-2][1:12]/
    #   sum(res.store[res.store[,ncol(res.store)-2] != 0 ,ncol(res.store)-2][1:12])
    # 
    # res.store[res.store[,ncol(res.store)-2] != 0 ,3][1:12]/
    #   sum(res.store[res.store[,ncol(res.store)-2] != 0 ,3][1:12])
    
    sel.m <- res$parent[res$parent$Sex == "male", ][
      order(res$parent[res$parent$Sex == "male","oc"], 
            decreasing = T)[1:n.sel.males], c("Indiv", "Sex", "oc")]
    Sire <- 
      as.integer(get.vector.occurence(sel.m, population.size = population.size))
    
    sel.f <- res$parent[res$parent$Sex == "female", ][
      order(res$parent[res$parent$Sex == "female","oc"], 
            decreasing = T)[1:n.sel.females], c("Indiv", "Sex", "oc")]
    Dam <- 
      as.integer(get.vector.occurence(sel.f, population.size = population.size))
    
    # sel.m
    # db.sel.m <- get.database(pop, id = sel.m[, "Indiv"])
    # bv.m <- get.bv(pop, database = db.sel.m, use.id = T)
    # 
    # db.sel.f <- get.database(pop, id = sel.f[, "Indiv"])
    # bv.f <- get.bv(pop, database = db.sel.f, use.id = T)
    # 
    # mean(get.bv(pop, gen = length(pop$breeding))[1,])
    # sum((sel.m[,"oc"]*bv.m)) +sum((sel.f[,"oc"]*bv.f))
    #res.i2 <- res
  }
  if(OGC.type == "no.OGC"){
    # no OGC means all animals that are selected based on bv
    # have the same probability to be drawn as a parent
    sex <- get.sex(pop, database = db)
    bvs <- get.bv(pop, database = db, use.id = T)
    
    ids.sel.m <- sort(as.integer(names(bvs[,sex == 1])[order(bvs[,sex == 1], decreasing = T)[1:n.sel.males]]))
    db.sel.m <- get.database(pop, id = ids.sel.m)
    
    ids.sel.f <- sort(as.integer(names(bvs[,sex == 2])[order(bvs[,sex == 2], decreasing = T)[1:n.sel.females]]))
    db.sel.f <- get.database(pop, id = ids.sel.f)
    
    db <- rbind(db.sel.m, db.sel.f)
    bvs <- get.bv(pop, database = db, use.id = T)[1,]
    
    sel.m <- data.frame(Indiv = ids.sel.m, oc = 1/length(ids.sel.m))
    Sire <- as.integer(get.vector.occurence(res = sel.m, 
                                 population.size = population.size))
    
    sel.f <- data.frame(Indiv = ids.sel.f, oc = 1/length(ids.sel.f))
    Dam <- as.integer(get.vector.occurence(res = sel.f, 
                                population.size = population.size))
  }
  
  ##############################################################################
  
  if(one.male.per.female){
    Sire <- get.vector.occurence(sel.m, population.size = n.sel.females)
    Dam <- sel.f[,"Indiv"]
    offspring.per.mating <- rep(population.size%/%n.sel.females, times = n.sel.females)
    pos <- sample(1:length(offspring.per.mating), replace = F, size = population.size%%n.sel.females)
    offspring.per.mating[pos] <- offspring.per.mating[pos]+1
    
  } else{
    offspring.per.mating <- 1
  }
  Mating <- data.frame(Sire, Dam, offspring.per.mating)
  Mating <- apply(Mating, MARGIN=2, as.integer)
  Mating
  
  m <- apply(Mating, MARGIN = 1, function(rowi) {
    list(matrix(
      rep(c(
        c(get.database(pop, id = rowi[1])[-3], 
          get.database(pop, id = rowi[2])[-3])
      ), rowi[3]), ncol = 6, byrow = T
    ))
  })
  m <- rlist::list.rbind(rlist::list.flatten(m))
  
  # this is randomly assigning the sex
  inaccuracy <- sample(c(-0.1,0.1),1) # this is added because if an uneven number 
  # is divided by 2, it gives 0.5 and then it is impossible to decide whether
  # there should be one more male or female
  m <- cbind(m, sample(c(
    rep(1, round((nrow(m)/2) + inaccuracy)),
    rep(0, round((nrow(m)/2) - inaccuracy))),
    nrow(m)))
  ##############################################################################
  
  pop <- breeding.diploid(pop,
                          fixed.breeding = m,
                          verbose = F,
                          delete.same.origin = T)
  
  k <- kinship.emp.fast(pop, gen = length(pop$breeding), ibd.obs = 1000, hbd.obs = population.size)
  vec.av.kin <- c(vec.av.kin, k[1])
  vec.av.inbreeding <- c(vec.av.inbreeding, k[2])
  vec.av.bv <- c(vec.av.bv, mean(get.bv(pop, gen = length(pop$breeding))[1,]))
  
  # need to do this otherwise the variance when having litters is weird
  pop2 <- breeding.diploid(pop, 
                           breeding.size = 500, 
                           selection.m.database = get.database(pop, gen = length(pop$breeding))[1,],
                           selection.f.database = get.database(pop, gen = length(pop$breeding))[2,], 
                           max.offspring = ceiling(c(500/(population.size/2), 500/(population.size/2))), 
                           max.mating.pair = 1, verbose = F)
  
  vec.gen.var <- c(vec.gen.var, var(get.bv(pop2, gen = length(pop2$breeding))[1,]))
  #get.effective.size(pop, gen = length(pop$breeding))
  cat(paste0("Finished generation ", i, "\n"))
}
################################################################################
cat("Finished simulation.")

# t1 <- Sys.time()
# k <- kinship.emp(pop, gen = length(pop$breeding)) #needs 1.5 minutes for 210 animals
# t2 <- Sys.time()
# t2-t1
# 
# kiemp <- k[upper.tri(k)]
# kexp <- kinship.exp(pop, gen = length(pop$breeding))
# kiex <- kexp[upper.tri(kexp)]
# 
# cor(kiemp, kiex)

if(OGC.type == "no.OGC"){
  valid.solution <- rep(NA, times = generations)
}

time.sim.end <- Sys.time()
runtime <- difftime(time.sim.end, time.sim.start, units = "min")[[1]]

out.table <- parameters
n.col <- length(c(runtime,vec.av.kin, vec.av.inbreeding, vec.av.bv, vec.gen.var, valid.solution))
out.table[rowi, (ncol(out.table)+1):(ncol(out.table)+n.col)] <- 
  c(runtime,vec.av.kin, vec.av.inbreeding, vec.av.bv, vec.gen.var, as.character(valid.solution))

colnames(out.table) <- c(colnames(parameters), 
                         "run time (minutes)",
                         paste0("mean kin gen ", 1:length(vec.av.kin)),
                         paste0("mean inbreeding gen ", 1:length(vec.av.inbreeding)),
                         paste0("mean bv gen ", 1:length(vec.av.bv)),
                         paste0("mean gen var gen ", 1:length(vec.gen.var)),
                         paste0("solution valid ", 1:length(valid.solution)))

time <- paste0("_",Sys.time())
time <- stringr::str_replace_all(time, pattern = '[-:]', "")
time <- stringr::str_replace_all(time, pattern = ' ', replacement = "_")

filename <- paste0(name_parameter_table,'_run',rowi,'_time-',time,'.txt')

data.table::fwrite(x = out.table, file = filename)

################################################################################
################################################################################
# this bit here is for merging the data so that I don't have to do it anymore
n_files_expected <- nrow(parameters)

ls_files <- list.files()
txt_files <- stringr::str_detect(ls_files, pattern = '\\.txt$')
txt_files <- ls_files[txt_files]

# numbers <- 0 is included to keep the remaining part executable even if no
# file was found
numbers <- 0
numbers <- stringr::str_extract_all(txt_files, pattern = "run[[:digit:]]+")
numbers <- stringr::str_extract_all(unlist(numbers), pattern = '[[:digit:]]+', simplify = FALSE)
numbers <- as.numeric(unlist(numbers))

numbers[duplicated(numbers)]

numbers <- unique(sort(numbers))

(1:n_files_expected)[!(1:n_files_expected %in% numbers)]

# the following should only become TRUE if all files are already written out
condition <- sum(numbers %in% 1:n_files_expected) == n_files_expected

if(condition){squ
  for(i in 1:n_files_expected){
    file <- txt_files[stringr::str_detect(txt_files, pattern = paste0('run',i,'_'))][1]
    if(is.na(file)){
      extracted <- stringr::str_extract_all(txt_files, pattern = paste0('run[[:digit:]]*_'))
      dummy_file <- str_extract(extracted[[2]], pattern = "[[:digit:]]+")
      file <- txt_files[stringr::str_detect(txt_files, pattern = paste0('run',dummy_file,'_'))][1]
    }
    IN <- data.table::fread(file)
    
    if(i == 1){out <- IN[i,]; next}
    
    out <- rbind(out, IN[i,], fill = TRUE)
  }
  # the same time as for the output above is used
  data.table::fwrite(out, file = paste0(name_parameter_table, '_all_time-',time,'.txt'))
  
  # here I am converting logical values to character
  ncol <- which(unlist(lapply(out, is.logical)))
  out <- as.data.frame(out)
  out[,ncol] <- unlist(lapply(out[,ncol], as.character))
  writexl::write_xlsx(out, paste0(name_parameter_table, '_all_time-',time,'.xlsx'))
}
################################################################################
