# Genetic algorithm optimization

reproduce_fixedpop <- function(popdata, modelparams) {

    # minimum to replace (based on preset replacement rate)
    minreplace <- ifelse(modelparams$replacerate == 1, modelparams$popsize - 1, modelparams$replacerate * modelparams$popsize)

    # indices of dead
    aliveandrepro <- popdata$alive & popdata$fitness >= 1
    dead <- which(!(aliveandrepro))

    # find number of dead to replace
    nreplace <- length(dead)

    # identify indices of eligible parent solutions
    parents <- which(aliveandrepro)

    # indices of solutions to be replaced
    children <- dead

    # if number to be replaced is greater than death count cull those with lowest fitness
    if (minreplace > nreplace) {
      lowest <- parents[order(popdata[parents]$fitness, decreasing = FALSE)][1:(minreplace - nreplace)]
      children <- c(children, lowest)
      nreplace <- minreplace
    }

    # SELECTION
    ## --------- Roulette Selection
    # choose new children using roulette selection where higher fitness == greater chance of being selected

    if (modelparams$selectiontype == "roulette") {
      roulette <- sample(parents,
        size = nreplace,
        prob = popdata[parents, fitness] / sum(popdata[parents, fitness]),
        replace = TRUE)
      selectedparents <- roulette
    }

    ## --------- Stochastic universal sampling (pg 50 in Sivanandam & Deepa 2008. Introduction To Genetic Algorithms)
    if (modelparams$selectiontype == "sus") {
      orderedparents <- parents[order(popdata[parents, fitness])]
      cumfit <- cumsum(popdata[orderedparents, fitness])
      last <- cumfit[length(cumfit)]
      interval <- last / nreplace
      pointers <- seq(from = runif(1, 0, interval), to = last, by = interval)
      susselect <- as.integer(findInterval(pointers, cumfit) + 1)
      selectedparents <- orderedparents[susselect]
    }

    if (modelparams$selectiontype == "rank") {
      ## --------- Rank selection
      ranktourney <- function(t) {
        tourney <- sample(parents, t, replace = TRUE)
        best <- which.max(popdata[tourney]$fitness)
        best <- tourney[best]
        return(best)
      }

      # defaults to tournament size of 2 (increase to make more proportionally representative)
      selectedparents <- replicate(nreplace, ranktourney(2))
    }

    # add clones to population
    popdata[children, (names(popdata)):= popdata[selectedparents,]]
    popdata[children, age := 0]

    # CROSSOVER

    # get traits for crossover/mutation, m_traits are time related whole numbers, e_traits are continuous
    e_traits <- modelparams$e_traits
    m_traits <- modelparams$m_traits
    all_traits <- c(e_traits, m_traits)

    # randomly identify crossed offspring
    crossed <- as.logical(rbinom(children, 1, modelparams$crossoverrate))
    # make number of crossed even
    if((sum(crossed) %% 2) > 0){
      crossed[sample(which(crossed == 0))] <- 1
    }

    # if any are crossed then cross them from 2 random parents in the parent pool
    # if (sum(crossed) > 1) {
      # get parent pairs for crossing using niching (probabilistically pair more holistically similar parent solutions)
      crossparents1 <- sample(selectedparents, sum(crossed)/2, replace = T)
      if(modelparams$niching){
        simmat <- as.matrix(proxy::simil(scale(popdata[unique(crossparents1), ..all_traits], center = F), diag = T, upper = T, method = "cosine"))
        crossparents2 <-
          sapply(crossparents1, FUN = function(x) sample(unique(crossparents1), size = 1,
                                                         prob = simmat[,unique(crossparents1) %in% x]/sum(simmat[,unique(crossparents1) %in% x])))
      } else {
        crossparents2 <- sapply(crossparents1, FUN = function(x) sample(selectedparents[!(selectedparents %in% x)], 1))
      }

     ### random blended crossover of energetic (continous) traits (random prop of parental solution contribution)
      # proportion of parent solution contribution
      prop <- runif((sum(crossed)/2),0,1)

      # child solutions to be crossed
      childrencross <- children[crossed]

      # divide crossed solutions into two equal pools
      cc1 <- childrencross[1:(length(childrencross)/2)]
      cc2 <- childrencross[(length(childrencross)/2 + 1):length(childrencross)]

      # energetic traits to blend
      e_blend <- e_traits

      # time/movement track traits to blend
      m_blend <- m_traits

      # cross and blend
      set(popdata, i = cc1, j = e_blend,
          value = popdata[crossparents1, ..e_blend] - (prop * (popdata[crossparents1, ..e_blend] - popdata[crossparents2, ..e_blend])))

      set(popdata, i = cc2, j = e_blend,
          value = popdata[crossparents2, ..e_blend] - (prop * (popdata[crossparents2, ..e_blend] - popdata[crossparents1, ..e_blend])))

      set(popdata, i = cc1, j = m_blend,
          value = floor(popdata[crossparents1, ..m_blend] - (prop * (popdata[crossparents1, ..m_blend] - popdata[crossparents2, ..m_blend]))))

      set(popdata, i = cc2, j = m_blend,
          value = floor(popdata[crossparents2, ..m_blend] - (prop * (popdata[crossparents2, ..m_blend] - popdata[crossparents1, ..m_blend]))))

    # MUTATION
    # identify which members of the population to mutate
     mutated <- as.logical(rbinom(nreplace, 1, modelparams$mutationrate))

    # mutate random portion of all traits and add mutated to population
    if (any(mutated)) {
      # mutate traits by blending with random prop of randomly initialized values
      propmut <- runif(sum(mutated),0,1)
      mutated <- children[mutated]

      set(popdata, i = mutated, j = e_traits,
          value = (1-propmut) * popdata[mutated, ..e_traits] + propmut * as.data.table(randomize_traits(modelparams,n = length(mutated)))[, ..e_traits])

      set(popdata, i = mutated, j = m_traits,
          value = floor((1-propmut) * popdata[mutated, ..m_traits] + propmut * as.data.table(randomize_traits(modelparams,n = length(mutated)))[, ..m_traits]))

    }
}

