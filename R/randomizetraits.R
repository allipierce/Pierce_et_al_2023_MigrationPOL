# Random initialization of IBM parameters

randomize_traits <- function(modelparams,
                             #world,
                             traits = NULL, n = NULL, popdata = NULL, ...){
  with(modelparams, {
    if(is.null(traits)) traits <- c(e_traits, m_traits)
    if(is.null(n)) n <- popsize
    if(!is.null(popdata)) {
      mapply(assign, e_traits, popdata[,..e_traits])
      mapply(assign, m_traits, popdata[,..m_traits])
    }

    # zoom factor - ratio of max body lengths of compared species
    # set to 1 for now to maintain relative simplicity
    z <- 1

    # somatic energy stream partition parameter
    kappa <- runif(n, kappa_min, kappa_max) #kappa from the deb model, proportion to somatic

    # minimum value for v
    minval <- 0.01

    # energy conductance parameter. Ratio between how much energy can be mobilized from reserves relative to how fast it is filled per time step. Smaller values indicate slower mobilization rates, energy it put into reserves faster than it is pulled out.
    v <- runif(n, minval, 0.99)

    # minimum lb defined by delta l at t0 = v/3
    l_b <- runif(n, v/3, 0.99)

    # g (maximum energy investment ratio) restricted such that potential energy of 1 unit of volume cannot exceed maximum energy density but is higher than the minimum energy density required to achieve birth (eb > l_b).
    # In the case of Lm = 1 and 100% growth efficiency that means g < 1/kappa and g > l_b/kappa
    # [EG] <= [Em]Lm
    # [EG]/[Em] <= [Em]/[Em]
    # g * kappa <= 1
    # [EG]/[Em] >= l_b
    # g * kappa >= l_b

    ming <- l_b/kappa
    maxg <- 1/kappa
    g <- runif(n, ming, maxg)

    #maturity threshold for puberty
    uph <- runif(n, (1-kappa) * l_b^3, (1-kappa))

    #a_b limits defined in Kooijman 2010 used for checking
      # max_ab <- (3 * l_b)/v * (1 + g)
      # min_ab <- (3 * l_b)/v
      #
      # max_ab <- (3 * l_b)/v * (1 + (g/l_b))
      # min_ab <- (3 * l_b)/v

    #calc scaled time at birth to calculate age at birth, round age to nearest whole number
    t_b <- calct_b_v(lb = l_b, g = g, eb = l_b)
    a_b <- t_b/(v/g)
    a_b <- as.integer(round(a_b))

    # Used for checking with a_b limits commented above
    # summary(max_ab - a_b)
    # summary(a_b - min_ab)
    # summary(a_b)

if(movement){
  # generate a non-breeding location duration, can't be longer than total time steps - shortest possible migration (4 time steps)
  nb_length <- as.integer(runif(n, 1, steps - 4))
  # generate duration of single season migration cannot be longer than half the difference between total time steps and non-breeding location duration.
  mig_time <- as.integer(runif(n, 1, (steps - nb_length)/2 ))

      # cost of transport per surface area per 1 unit of distance rel to volume specific maint rate
      m <- runif(n, 0.01, 1)

    } else {
      mig_time <- 0
      nb_length <- 0
      m <- 0
    }
    gc()
    rtraits <- mget(traits)
    return(rtraits)
  })
}

initialize_pop <- function(popdata, modelparams
                           ){
  # set start and end locations, set solution age to 0 (solution age for internal tracking/tuning purposes)
  popdata[, `:=`(startloc = modelparams$startloc,
                 endloc = modelparams$endloc,
                 age = 0)
  ]

  # sets initial states
  reset_pop(popdata, modelparams)
}

# function to reset starting states
reset_pop <- function(popdata, modelparams){
  popdata[, `:=`(
                 startloc = modelparams$startloc,
                 endloc = modelparams$endloc,
                 loc = modelparams$startloc,
                 # starting reserves assume maternal effect (reserves of mother at laying same as embryo at hatch) and that reserves of mother are equal to energy at birth location at time of birth
                 reserves = lookupf_bylat((modelparams$steps - a_b), modelparams$startloc, modelparams),
                 distance = 0,
                 nb_start = 0,
                 repro_reserves = 0,
                 alive = TRUE,
                 t_mat = as.numeric(NA),
                 # km value based on Lm = 1 for simplicity in limiting g (would scale with z if implemented)
                 km = v/g,
                 mig_speed = if(modelparams$movement){
                   abs(modelparams$startloc - modelparams$endloc)/mig_time
                 }else{0},
                 ubh = (1 - kappa) * l_b^3,
                 uh = (1 - kappa) * l_b^3,
                 l = l_b,
                 l_p = as.numeric(NA),
                 fitness = 0,
                 eb = as.numeric(NA),
                 eggcost = as.numeric(NA)
                 )]
}

