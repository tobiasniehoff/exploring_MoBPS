# returns the object needed for fixed mating
# population is the population object used in mobps
# id_m is vector with ids of males like c(32,12,4)
# n.offspring is number of offspring
# this function assumes that within sex, the contributions per animal are equal
# (e.g., every male has the same number of offspring)
#
# The function uses the A matrix based on pedigree info only to get 
# relationships.
#
# in my tests, the function optiSel::matings() has often not found an acceptable 
# solution. I find this especially the case when many males and females are selected
# But I still want consistent output. So when mate allocation of optiSel fails,
# I assign mates randomly
# I assume that equally as many females as males are produced
#
# the code will first try to do mate allocation with solver lpsymphony 
# if lpsymphony is not available, it will try to do it with Rsymphony,
# If Rsymphony is not available, it wil ltry a default solver
# if this is not working, then it will try the default solver in optiSel
# this takes long and does not always provide a solution
# if this also fails, then mates are allocated so that fullsib and halfsib
# matings are avoided
# to directly allocate mates with avoiding halfisb and fullsib matings without
# first trying a solver, set avoid_sibmating to TRUE
#
#
min.coanc.mating <- function(population, n.offspring, id_m, id_f,
                             avoid_sibmating = F, depth.pedigree = 7, 
                             time.limit.sec = 3600){
  library(MoBPS)
  # # ids of males
  # id_m <- get.id(population, 
  #                cohorts = get.cohorts(population)[c(T,F)][gen], 
  #                use.id = T)
  # 
  # # ids of females
  # id_f <- get.id(population, 
  #                cohorts = get.cohorts(population)[c(F,T)][gen], 
  #                use.id = T)
  
  Indiv <- c(id_m, id_f)
  Sex <- c(rep("male" ,length(id_m)), rep("female" ,length(id_f)))
  
  n.off.male <- rep(n.offspring%/%length(id_m), length(id_m))
  n.off.male[sample(1:length(id_m), n.offspring%%length(id_m))] <- 
    n.off.male[sample(1:length(id_m), n.offspring%%length(id_m))] + 1
  
  n.off.female <- rep(n.offspring%/%length(id_f), length(id_f))
  n.off.female[sample(1:length(id_f), n.offspring%%length(id_f))] <- 
    n.off.female[sample(1:length(id_f), n.offspring%%length(id_f))] + 1
  
  n <- c(n.off.male, n.off.female)
  
  candi <- data.frame(Indiv, Sex, n)
  if(!avoid_sibmating){
    
    Amat <- kinship.exp(population, 
                        database = get.database(population, id = c(id_m, id_f)),
                        depth.pedigree = depth.pedigree)
    
    if("lpsymphony" %in% installed.packages()){
      Mating <- R.utils::withTimeout(
        optiSel::matings(phen = candi,
                         Kin = Amat,
                         solver = lpsymphony::lpsymphony_solve_LP,
                         ub.n = NA, # this says there is no constraint for how often an individual can be 
                         # mated to the same other individual
                         max = FALSE # means objective function is minimize -> low inbreeding
        ), timeout = time.limit.sec, onTimeout = "warning"
      )
      # https://stackoverflow.com/questions/34346619/how-to-stop-a-function-in-r-that-is-taking-too-long-and-give-it-an-alternative
    } else if("Rsymphony" %in% installed.packages()){
      Mating <- optiSel::matings(phen = candi,
                                 Kin = Amat,
                                 solver = Rsymphony::Rsymphony_solve_LP,
                                 ub.n = NA, # this says there is no constraint for how often an individual can be 
                                 # mated to the same other individual
                                 max = FALSE # means objective function is minimize -> low inbreeding
      )
    } else{
      Mating <- optiSel::matings(phen = candi,
                                 Kin = Amat,
                                 solver = "default",
                                 ub.n = NA, # this says there is no constraint for how often an individual can be 
                                 # mated to the same other individual
                                 max = FALSE # means objective function is minimize -> low inbreeding
      )
      
      cat("Consider installing package lpsymphony or Rsymphony for mate allocation and faster run time. \n")
    }
  } 
  # it is possible that no ideal solution is found for minimum coancestry mating
  # in that case, I do random mating
  # thought: allocating matings fails sometimes with optiSel. Why does it work
  # in Jibrila then/ what is the advantage of the software they are using?
  if(!avoid_sibmating){
    avoid_sibmating <- F
    if(!exists("Mating")){
      avoid_sibmating <- T
    } else if(!avoid_sibmating){
      if(is.na(Mating[1,1]) | nrow(Mating) == 0 | sum(Mating$n) != n.offspring){
        avoid_sibmating <- T
      }
    }
  }
  
  if(avoid_sibmating){
    
    # if the mate allocation of optiSel has failed, I need to do it myself
    # the solution is to allocate mates randomly. But an improvement of that 
    # is to avoid mating fullsibs and halfsibs
    if(TRUE){
      p <- MoBPS::breeding.diploid(
        population, 
        breeding.size = n.offspring,
        selection.m.database = get.database(population, id = id_m),
        selection.f.database = get.database(population, id = id_f),
        max.offspring = c((n.offspring/length(id_m)), 
                          (n.offspring/length(id_f))),
        avoid.mating.fullsib = T,
        avoid.mating.halfsib = T,
        verbose = F
      )
      Mating <- get.pedigree(p, gen = length(p$breeding), raw = F, id = T)[,1:3]
      Mating <- Mating[,c(2,3,1)]
      Mating[,3] <- 1
      Mating <- as.data.frame(Mating)
      colnames(Mating) <- c("Sire", "Dam" ,"n")
      
      cat("Minimum coancestry mating allocation failed.\n
          Doing random mate allocation but avoiding fullsib and halfsib mating.\n")
      
    } else{
      resample <- function(x, ...) {x[sample.int(length(x), ...)]} # idea from ?sample
      Sire <- resample(c(rep(id_m, n.offspring%/%length(id_m)),
                         sample(id_m, n.offspring%%length(id_m))))
      Dam <- resample(c(rep(id_f, n.offspring%/%length(id_f)),
                        sample(id_f, n.offspring%%length(id_f))))
      n <- rep(1, length(Dam))
      
      Mating <- data.frame(Sire, Dam, n)
      
      cat("Minimum coancestry mating allocation failed. Doing random allocation.\n")
    }
  } else{
    Mating <- apply(Mating, MARGIN=2, as.integer)
  }
  
  # with m I am trying to create the fixed.breeding object for Mobps
  m <- apply(Mating, MARGIN = 1, function(rowi) {
    list(matrix(
      rep(c(
        c(get.database(population, id = rowi[1])[-3], 
          get.database(population, id = rowi[2])[-3])
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
  #m
  
  return(m)
}
# population.db <- get.database(population, gen = length(population$breeding))
# id_m <- sample(get.id(population, database = population.db[population.db[,2]==1,]),5)
# id_f <- sample(get.id(population, database = population.db[population.db[,2]==2,]),20)
# 
# mate.plan <- min.coanc.mating(population = population,
#                               n.offspring = 200,
#                               id_m = id_m, id_f = id_f, 
#                               avoid_sibmating = FALSE, 
#                               depth.pedigree = 7)
# 
# population <- breeding.diploid(population,
#                         fixed.breeding = mate.plan,
#                         verbose = F,
#                         delete.same.origin = T)


# this function returns either a pedigree (return_database = F) or a database
# with all the ancestors tracing back to the desired number of generations
# most of the code is the source code from MoBPS::kinship.exp()
# id=TRUE is a bit faster
get.direct.ancestors <- function(population, 
                                 database = NULL, 
                                 gen = NULL, 
                                 cohorts = NULL,
                                 depth.pedigree=7, 
                                 return_database = F,
                                 elements = NULL,
                                 id = F,
                                 storage.save=1.5,
                                 verbose = T){
  if(depth.pedigree == 0){
    cat("Since depth.pedigree is set to 0, ancestors cannot be looked up.\n")
    return(NULL)
  }
  #if(depth.pedigree >= length(population$breeding)){depth.pedigree <- length(population$breeding)-1}
  if(depth.pedigree == Inf) {depth.pedigree <- length(population$breeding)}
  
  database <- get.database(population, gen=gen, database=database, cohorts=cohorts)
  
  n.animals <- sum(diff(t(database[,3:4, drop=FALSE]))+1)
  
  if(length(elements)>0){
    
    if(max(elements)> n.animals){
      stop("kinship emp number of individuals does not match!")
    }
    
    cumorder <- cumsum(c(1,diff(t(database[,3:4, drop=FALSE]))+1))
    activ_database <- rep(FALSE, nrow(database))
    elements_new <- elements
    for(index in 1:nrow(database)){
      if(sum(intersect(elements, cumorder[index]:(cumorder[index+1]-1)))>0){
        activ_database[index] <- TRUE
      } else{
        elements_new[elements>cumorder[index]] <- elements_new[elements>cumorder[index]] - database[index,4] + database[index,3] - 1
      }
    }
    
    elements <- elements_new
    database <- database[which(activ_database),,drop=FALSE]
    
    if(TRUE){
      cumorder <- cumsum(c(1,diff(t(database[,3:4, drop=FALSE]))+1))
      elements_new <- elements
      for(index in 1:nrow(database)){
        remain <- intersect(elements, cumorder[index]:(cumorder[index+1]-1)) - cumorder[index] +1
        
        if(max(remain) < (database[index,4] - database[index,3]+1)){
          elements_new[elements>(cumorder[index+1]-1)] <- elements_new[elements>(cumorder[index+1]-1)] + max(remain) - (database[index,4] - database[index,3] +1 )
          database[index,4] <- max(remain) + database[index,3] - 1
          
        }
      }
      elements <- elements_new
    }
    elements <- elements_new
    
  } else{
    elements <- 1:sum(database[,4]-database[,3]+1)
  }
  
  if(depth.pedigree==Inf){
    pedigree.database <- get.database(population, gen=1:max(database[,1]))
  } else{
    new.pedigree.database <- pedigree.database <- database
    remaining.depth <- depth.pedigree
    while(remaining.depth>0){
      parents <- get.pedigree(population, database = new.pedigree.database, raw=TRUE)
      m_parents <- rbind(parents[parents[,5]==1,4:6], parents[parents[,8]==1,7:9])
      f_parents <- rbind(parents[parents[,5]==2,4:6], parents[parents[,8]==2,7:9])
      if(nrow(m_parents)>0){
        m_gen <- unique(m_parents[,1])
        m_data <- cbind(m_gen, 1, 0,0)
        nincluded <- numeric(length(m_gen))
        for(index in 1:length(m_gen)){
          m_data[index,3] <- min(m_parents[m_parents[,1]==m_gen[index],3])
          m_data[index,4] <- max(m_parents[m_parents[,1]==m_gen[index],3])
          nincluded[index] <- length(unique(m_parents[m_parents[,1]==m_gen[index],3]))
        }
        
        for(index in length(m_gen):1){
          if(nincluded[index] < (m_data[index,4]-m_data[index,3]+1)/storage.save){
            m_data <- m_data[-index,]
            activ_p <- unique(m_parents[m_parents[,1]==m_gen[index],3])
            m_data <- rbind(m_data, cbind(m_gen[index], 1, activ_p, activ_p))
          }
        }
        
      } else{
        m_data <- NULL
      }
      if(nrow(f_parents)>0){
        f_gen <- unique(f_parents[,1])
        f_data <- cbind(f_gen, 2, 0,0)
        nincluded <- numeric(length(f_gen))
        for(index in 1:length(f_gen)){
          f_data[index,3] <- min(f_parents[f_parents[,1]==f_gen[index],3])
          f_data[index,4] <- max(f_parents[f_parents[,1]==f_gen[index],3])
          nincluded[index] <- length(unique(f_parents[f_parents[,1]==f_gen[index],3]))
        }
        
        for(index in length(f_gen):1){
          if(nincluded[index] < (f_data[index,4]-f_data[index,3]+1)/storage.save){
            f_data <- f_data[-index,]
            activ_p <- unique(f_parents[f_parents[,1]==f_gen[index],3])
            f_data <- rbind(f_data, cbind(f_gen[index], 2, activ_p, activ_p))
          }
        }
        
      } else{
        f_data <- NULL
      }
      
      new.pedigree.database <- get.database(population, database=rbind(m_data,f_data))
      new.pedigree.database <- unique(new.pedigree.database)
      remaining.depth <- remaining.depth - 1
      pedigree.database <- rbind(new.pedigree.database, pedigree.database)
    }
    
    pedigree.database <- get.database(population, database = pedigree.database)
  }
  
  ids_database <- get.id(population, database = database)
  ids_database_unique <- unique(ids_database)
  ids_pedigree <- sort(unique(get.id(population, database = pedigree.database)))
  
  ids_pedigree_first <- max(get.id(population, 
                                   database = pedigree.database[pedigree.database[1,1]==pedigree.database[,1],,drop=FALSE]))
  
  n.animals <- length(ids_database_unique)
  n.total <- length(ids_pedigree)
  
  if(verbose) cat("Derive pedigree for ", n.animals, " individuals based on ", n.total, " individuals.\n")
  
  if(return_database == F){
    ped <- get.pedigree(population = population, 
                        database = pedigree.database, id = id)
    return(ped)
  }
  if(return_database == T){
    ped <- get.pedigree(population = population, 
                        database = pedigree.database, id = T)
    id.focal.inds <- get.id(population = population, database = database, use.id = T)
    # if any individual in the focal group is parent to another individual, it needs
    # to be kept. later, all the focal inds are thrown out of the pedigree
    id.focal.inds <- id.focal.inds[!(id.focal.inds %in% c(ped[,2], ped[,3]))]
    
    ped <- ped[!(ped[,1] %in% id.focal.inds),]
    db <- get.database(population = population, id = as.integer(ped[,1]))
    
    return(db)
  }
}

# population <- creating.diploid(nsnp = 1000, chr.nr = 30, nindi = 100,
#                         n.additive = 1000, base.bv = 100)
# # some selection
# for(i in 1:8){
#   population <- breeding.diploid(population,
#                    selection.size = c(5,20),
#                    breeding.size = 100,
#                    selection.criteria = c("bv", "bv"),
#                    selection.m.cohorts = get.cohorts(population)[length(get.cohorts(population))-1],
#                    selection.f.cohorts = get.cohorts(population)[length(get.cohorts(population))])
# }
# 
# get.direct.ancestors(population=population, 
#                      cohorts = "Cohort_6_F", 
#                      depth.pedigree = 2)