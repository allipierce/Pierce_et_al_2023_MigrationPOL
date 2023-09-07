## DEB functions

#assimilation function
#amount of scale energy assimilated is a function of surface area (l^2)
assimilate <- function(l, f) {
  f * l^2
}

# mobilization function
# amount mobilized into the body depends on energy reserves, current length and investment ratio
mobilize <- function(e, l, g, l_t = 0, l_move = 0) {
  e * l^2 * ((g + l + l_t + l_move)/(g + e))
}

# somatic maintenance function
# amount allocated to maintenance of somatic processes which depends on volume. l_t = 0 if environmental temp is constant
maintain_soma <- function(kappa, l, l_t = 0, l_move = 0) {
  kappa * l^2 * (l + l_t + l_move)
}

# maturity maintenance function
# energy allocated to maintaining reproductive status calculated as a fraction of scaled cost of maturity (fraction is cost to maturity for each step, default is 1 - maturity maintenance equal to somatice maintence)
maintain_maturity <- function(k = 1, uh) {
  k*uh
 }

# reproduction function
# energy allocated to maturity if immature or reproductive buffer if mature
reproduce <- function(kappa, e, l, g, l_t = 0, l_move = 0, k = 1, uh = uh){
  (1-kappa) * (e * l^2) * ((g + l + l_t + l_move)/(g + e)) - (k * uh)
}

# growth cost function, when e >= l this is equal to (kappa * mobilize) - maintain_soma
# energy allocated to growth based on current scaled length, reserves, and energy investment ratio
growcost <- function(kappa, l, e, g, l_t = 0, l_move = 0) {
  kappa * l^2 * ((e - l - l_t - l_move)/(1 + (e/g)))
}

#convert energy to length (grow)
grow <- function(p_g, l, g, kappa, km){
  p_g * km / (3 * l^2 * kappa)
}

expendenergy <- function(popdata, k, modelparams){

  # get available energy for that location
  c_loc <- popdata[,loc]

  # get f for current location
  f <- lookupf_bylat(k, c_loc, modelparams)

  # unpack DEB compound parameters
  l <- popdata[,l]
  l_p <- popdata[,l_p]
  e <- popdata[,reserves]
  ur <- popdata[,repro_reserves]
  g <- popdata[,g]
  kappa <- popdata[,kappa]
  m <- popdata[,m]
  km <- popdata[,km]
  uh <- popdata[,uh]
  uph <- popdata[,uph]
  z <- popdata[,z]
  distance <- popdata[,distance]/(z)
  t_mat <- popdata[,t_mat]
  a_b <- popdata[,a_b]
  eb <- popdata[,eb]

  # calc non-feeding movement costs (locomotion costs scale with surface area but ultimately scale with surface area due to distance scaling with body size)
  l_move <- m * distance

  # calculate what proportion of avail energy is assimilated based on scaled body size
  p_a <- assimilate(l = l, f = f)

  # mobilize from reserves based current scaled length and investment ratio
  # this is calculated in parts at p_s and p_g (p_c = p_s + p_g)
  # equal to p_s + p_g + p_r
  p_c <-  mobilize(e = e, l = l, g = g, l_move = l_move)

  # calc growth costs to see if starving and to increment length
  p_g <- growcost(kappa = kappa, l = l, e = e, g = g, l_move = l_move)
  #p_g <- (kappa * p_c) - p_s
  #p_g <- ifelse(p_g < 0, 0, p_g)
  #p_g <- ifelse(l == 1, 0, p_g)

  # calc repro costs to increment maturity or repro buffer and adjust for starvation
  p_repro <- reproduce(kappa = kappa, l = l, e = e, g = g, uh = uh, l_move = l_move)

  # if starving, subtract deficit in p_g from p_repro
  p_repro <- ifelse(p_g < 0, p_repro + p_g, p_repro)
  p_g <- ifelse(p_g < 0, 0, p_g)

  # if p_repro negative (not enough assimilated to cover maint and growth) no maturation (maturity can't be reversed but can be delayed)
  duh <- ifelse(p_repro < 0, 0, (p_repro * km))
  duh <- ifelse(uh >= uph, 0, duh)

  # calc new maturity level
  newuh <- uh + duh

  # if maturity meets puberty threshold then fix maturity to puberty threshold level
  newuh <- ifelse(newuh >= uph, uph, newuh)

  # if repro negative adjust to 0 to avoid depleting from repro-buffer
  p_repro_adj <- ifelse(p_repro < 0, 0, p_repro)

  # if mature, divert reproductive energy from maturation to reproduction
  # scale by km only not km*g*l^3 because we don't want repro energy buffer to be in terms of energy density [Er]/[Em] but in energy scaled by length only so that it can be later divided by estimated egg cost at the end

  dur <- ifelse(uh < uph, 0, (km * p_repro_adj))
  # don't tally up repro-buffer after lay date
  ur <- ifelse(k <= ((modelparams$step) - a_b), dur + ur, ur)

  # calc change in reserves and convert to units of energy density
  # if growth is zero don't adjust delta e for change in length (no shrinking)
  de <- g*km/l^3 * (p_a - p_c - (e/(g*kappa) * p_g))

  # adjust reserves - if repro reserves were insufficient to meet deficit, set e to -1 to flag as dead
  # currently not set to allow for use of repro-reserves
  e <- ifelse(ur < 0, -1, e + de)
  e <- ifelse(e < 0, -1, e)

  #if not starving - grow (if p_g is 0 then will not grow)
  newl <- l + grow(p_g, l, g, kappa, km)

  # if matured this iteration record length at maturity
  l_mat <- ifelse(newuh >= uph & is.na(t_mat), newl, l_p)

  # if matured this iteration record time of maturity
  newt_mat <- ifelse(newuh >= uph & is.na(t_mat), k, t_mat)

  # if current step is time of egg formation for offspring to be a_b at time = 1, record reserve for later checks of offspring viability
  neweb <- ifelse(k == ((modelparams$step) - a_b), e, eb)

  # record new states
  popdata[, `:=` (reserves = e,
                  eb = neweb,
                  l = newl,
                  repro_reserves = ur,
                  uh = newuh,
                  t_mat = newt_mat,
                  l_p = l_mat
                  )]

  # mark dead and zero out their repro buffer
  set(popdata, i = which(e < 0), j = "alive", value = F)
  set(popdata, i = which(!popdata$alive), j = "repro_reserves", value = 0)

  # cleanup unpacked DEB compound parameters
  rm(c_loc,f,l,l_p,e,ur,g,kappa,m,km,uh,uph,z,distance,t_mat,a_b,eb)
}
