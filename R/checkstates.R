## Functions to check viability of solutions

## Check if initial age and length at birth are viable
check_initviability <- function(popdata, modelparams){
  #mark dead if age at birth is greater than 1 year
  popdata[a_b > modelparams$steps, alive := FALSE]
  #mark dead if initial reserves less than length at birth
  popdata[reserves < l_b, alive := FALSE]
}

## Check if living and reproductive
check_livingandrepro <- function(popdata){
  #mark dead if reserves depleted
  popdata[reserves < 0, alive := FALSE]
  #mark dead if reproductive buffer empty
  popdata[repro_reserves <= 0, alive := FALSE]
  #mark dead if sexual maturity not reached
  popdata[uh < uph, alive := FALSE]
  #set e at time of birth to 0 if initial a_b was not possible
  popdata[is.na(eb), eb := 0]
  #mark dead if insufficient reserves for offspring to be born (growth ceases during development if embryo reserves at birth less than length at birth)
  popdata[eb < l_b, alive := FALSE]
}

## Check if completed migration in time
check_fullmig <- function(popdata, modelparams){
  # get current location
  c_loc <- popdata$loc
  # check if current location is close enough to starting location (within 1 movement)
  fullmig <- abs(c_loc - popdata$startloc) < popdata$mig_speed
  # mark dead if not close enough to breeding location
  popdata[!fullmig, alive := FALSE]
}

check_birthtiming <-  function(popdata, modelparams){
  # get alive and reproductive
  aliveandrepro <- which(popdata$alive == T & popdata$repro_reserves > 0)

  # calc age at birth of offspring given reserves at reproduction
  t_b <- calct_b_v(g = popdata[aliveandrepro]$g,
                   lb = popdata[aliveandrepro]$l_b,
                   eb = popdata[aliveandrepro]$eb)

  # get new ab based on eb at time of conception
  newa_b <- as.integer(round(t_b/popdata[aliveandrepro]$km))

  # check if new a_b is valid (> 0) and is smaller that number of steps, and is possible given maturation time, migration time, and non-breeding duration
  timefordev <- (newa_b >= 0) & (newa_b < modelparams$steps) & ((modelparams$steps - (popdata[aliveandrepro]$t_mat + (2 * popdata[aliveandrepro]$mig_time) + popdata[aliveandrepro]$nb_length)) >= newa_b)

  # if not enough time, kill solution
  popdata[aliveandrepro[!timefordev], alive := FALSE]

  # save new a_b to valid solutions
  popdata[aliveandrepro[timefordev], a_b := newa_b[timefordev]]
}
