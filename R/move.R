# Movement model function

fixed_move <- function(popdata, x, modelparams) {
  # get age of maturity
  t_mat <- popdata[,t_mat]
  # get current location
  c_loc <- popdata[,loc]
  #c_loc <- popdata[,sapply(loc, FUN = function(l) l[x])]
  # get end point
  endloc <- popdata[,endloc]
  # get start point
  startloc <- popdata[,startloc]
  # get duration of non-breeding period
  nb_length <- popdata[,nb_length]
  # get start of non-breeding period
  nb_start <- popdata[,nb_start]
  # get migrations speed
  mig_speed <- popdata[,mig_speed]

  # assess if it's time to move - mature and non-breeding period has not started, or current time is greater than end of non-breeding period and current location is farther than 1 movement from end location
  move <- (!is.na(t_mat)) & (nb_start == 0 | (x > (nb_start + nb_length) & abs(c_loc - startloc) >= mig_speed))

  # if moving, step length is same as migration speed, otherwise it is 0
  steplength <- ifelse(move, mig_speed, 0)

  # record new location - adjust direction of movement for hemisphere and direction of migration
  newloc <- ifelse(startloc > endloc,
                       ifelse(nb_start == 0, c_loc - steplength, c_loc + steplength),
                       ifelse(nb_start == 0, c_loc + steplength, c_loc - steplength))

  # if arrived at non-breeding location, record time non-breeding period started
  newnb_start <- ifelse(nb_start == 0 & abs(newloc - endloc) < mig_speed,
                     x,
                     nb_start)

  # record distance traveled
  newdistance <- steplength

  # update parameters for energy calculations
  popdata[, `:=` (loc = newloc,
                   nb_start = newnb_start,
                   distance = newdistance
                  )]

  # record the track
  popdata[,paste0("step",x) := newloc]

}



