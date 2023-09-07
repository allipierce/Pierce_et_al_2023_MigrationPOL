# Function to generate f based on location in a seasonal environment where resources vary sinusoidal by space and oscillate with time
# means can also be varied linearly over space with meanrange argument
# meanf > 0.5 will raise lower limit of f to keep f max at 1
# meanf < 0.5 will lower upper limit of f to keep f min > 0

lookupf_bylat <- function(t,loc, modelparams){
  with(modelparams,{
    mean_f <- sort(meanrange)
    mean_f <- mean_f[2] - ((mean_f[2]-mean_f[1]) * abs(loc)/(worldsize/2))
    amp <- loc/(worldsize/2) * (1 - mean_f)
    f <- amp * cos(2 * pi * t/steps) + mean_f
    return(f)
  })
}

