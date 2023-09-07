# GAEMM Model
# optimize metabolic (M) and movement (M) IBM parameters with a genetic algorithm (GA)
# data.table used to allow for more efficient memory management (can rewrite values by reference)

GAEMM <- function(modelparams, verbose = TRUE, lastonly = FALSE, tuning = FALSE, niching = TRUE, initdata = NULL) {
  # restrict data.table to single threaded internal processes (to reduce load during parallelization)
  data.table::setDTthreads(threads = 1L)
  # generate initial parameter values
  pdata <- as.data.table(randomize_traits(modelparams))
  # set initial states
  initialize_pop(pdata, modelparams)
  if(!is.null(initdata)) {
      traits <- c(modelparams$e_traits, modelparams$m_traits)
      set(pdata, i = 1:nrow(initdata), j = traits,
          value = initdata[, ..traits])
      reset_pop(pdata, modelparams)
    }

  # generation loop
  for (i in 1:modelparams$generations) {

    if(i > 1){
      # if beyond first generation use maternal reserves at birth from previous generation to initialize reserves
      preveb <- pdata[,eb]
      reset_pop(pdata, modelparams)
      pdata[, reserves := preveb]
    }

    check_initviability(pdata, modelparams)

    if(all(!pdata$alive)) {
      print(paste("population extinct || ", "generation: ", i, "step: ", k))
      break
    }

    # IBM-DEB movement model loop]
    for (k in 1:modelparams$steps) {
      expendenergy(pdata, k, modelparams)

      # move
      if(modelparams$movement){
      fixed_move(pdata, k, modelparams = modelparams)
      } else {
        # if no movement record current loc as starting location
        pdata[,c(paste0("step",k)) := modelparams$startloc]
      }

      # If all are dead then quit model run
      if(all(!pdata$alive)) {
        print(paste("population extinct || ", "generation: ", i, "step: ", k))
        break
      }
    }

    # mark any solutions with depleted reserves or are non-reproductive as dead
    check_livingandrepro(pdata)

    # if able to move then check that solutions made their migration
    if(modelparams$movement) check_fullmig(pdata, modelparams)

    # mark solutions as dead if offspring inviable
    check_birthtiming(pdata, modelparams)

    # calculate egg costs and fitness, unfit solutions marked as dead
    pdata[alive == T , eggcost := Vectorize(get_u0)(eb = eb, lb = l_b, g = g)]
    pdata[alive == T, fitness := as.integer(round(repro_reserves / eggcost))]
    pdata[fitness < 1, alive := FALSE]
    pdata[fitness < 1, fitness := 0]

    # increment internal solution age tracker
    pdata[, age := age + 1]

    # GA Model - Reproduction

    # if this is not the last loop then proceed with reproduction
    if(verbose){
      print(paste0("Generation: ", i))
    }

    # return parameters and model tuning stats if tuning run
    if(isTRUE(tuning)){
      if(i == 1) results <- data.frame()
      parmnames <- c(modelparams$e_traits, modelparams$m_traits)
      results[i,parmnames] <- pdata[alive & age > 0][which.max(fitness), ..parmnames]
      results$meanfitness[i] <- pdata[pdata$alive & age > 0, mean(fitness)]
      results$maxfitness[i] <- pdata[pdata$alive & age > 0, max(fitness)]
      results$varfitness[i] <- pdata[pdata$alive & age > 0, var(fitness)]
      results$nliving[i] <- pdata[pdata$alive & age > 0, .N]
      results$nunique[i] <- sum(!duplicated(data.frame(pdata[pdata$alive & age > 0])))
      results$generation[i] <- i
    } else {
      if(isTRUE(lastonly)){
        results <- as.data.frame(pdata)
      } else {
        pdata[, generation := i]
        results <- rbind(results, pdata)
      }
    }

    # break outer generation loop if extinct
    if(all(!pdata$alive)){
      print(paste("population extinct || ", "generation: ", i))
      break
    }

    if (i < modelparams$generations) {

      # reproduce (dead will be replaced in this step)
      reproduce_fixedpop(pdata, modelparams)
      # diagnostic readout of fitness
      if(verbose){
        print(paste0(
          "Mean Fitness: ", mean(pdata[pdata$alive & fitness >= 0]$fitness), " | Median Fitness: ", median(pdata[pdata$alive & fitness >= 0]$fitness)," | Max Fitness: ",
          max(pdata[pdata$alive & fitness > 0]$fitness), " | Fitness variance: ",var(pdata[pdata$alive & fitness > 0]$fitness)
        ))
        print(paste0("Popsize: ", nrow(pdata[pdata$alive])))
      }
    }
  }
  # diagnostic readout of fitness
  if(verbose){
    print(paste0(
      "Mean Fitness: ", mean(pdata[pdata$alive & fitness >= 0]$fitness), " | Median Fitness: ", median(pdata[pdata$alive & fitness >= 0]$fitness)," | Max Fitness: ",
      max(pdata[pdata$alive & fitness > 0]$fitness), " | Fitness variance: ",var(pdata[pdata$alive & fitness > 0]$fitness)
    ))
    print(paste0("Popsize: ", nrow(pdata[pdata$alive])))
  }
    return(results)
}

# function to run model with fixed parameters
GAEMMfixedparams <- function(pdata, modelparams){
  # preveb <- pdata$eb
  pdata <- as.data.table(pdata)
  # initialize_pop(pdata, modelparams)
  # pdata[,reserves := preveb]

  preveb <- pdata[,eb]
  reset_pop(pdata, modelparams)
  pdata[, reserves := preveb]
  check_initviability(pdata, modelparams)

  for (k in 1:modelparams$steps) {
    expendenergy(pdata, k, modelparams)
  }

  # mark any solutions with depleted reserves or are non-reproductive as dead
  check_livingandrepro(pdata)

  # mark solutions as dead if offspring inviable
  check_birthtiming(pdata, modelparams)

  # calculate egg costs and fitness, unfit solutions marked as dead
  pdata[alive == T , eggcost := Vectorize(get_u0)(eb = eb, lb = l_b, g = g)]
  pdata[alive == T, fitness := as.integer(round(repro_reserves / eggcost))]
  pdata[fitness < 1, alive := FALSE]
  pdata[fitness < 1, fitness := 0]
  return(pdata)
}
