# Code adapted from debtool matlab functions written by Bas Kooijman available at https://github.com/add-my-pet/DEBtool_M/tree/master/animal

get_beta0 = function(x0,x1){

  x03 = x0 ^ (1/ 3)
  x13 = x1 ^ (1/ 3)
  a3 = sqrt(3)
  f1 = - 3 * x13 + a3 * atan((1 + 2 * x13)/ a3) - log(as.complex(x13 - 1)) + log(as.complex(1 + x13 + x13 ^ 2))/ 2
  f0 = - 3 * x03 + a3 * atan((1 + 2 * x03)/ a3) - log(as.complex(x03 - 1)) + log(as.complex(1 + x03 + x03 ^ 2))/ 2
  Re(f1 - f0)
}

get_u0 <- function(eb = 1, lb, g){
  xb = g/ (eb + g)
  uE0 = (3 * g/ (3 * g * xb^(1/ 3)/ lb - get_beta0(0, xb)))^3
  return(uE0)
}

# calc tau_b to get a_b, adapted from debtool code. Uses integrate instead of pracma::quad since it can be marginally faster and handles infinite intervals, can be less accurate in some situations but testing didn't indicate it was in this case
calct_b <- function(lb, g, eb = 1){
  xb <- g/(eb + g)
  if(is.na(xb)) return(-1)
  ab <- 3 * g * xb^(1/3)/ lb
  3 * integrate(dget_tb, 1e-15, xb, ab = ab, xb = xb)$value
}

# integrand - adapted from debtool code
dget_tb <- function(x, ab, xb){
  # called by get_tb
  x ^ (-2/ 3) / (1 - x) / (ab - get_beta0(x, xb))
}

# vectorizes integrand function
calct_b_v <- Vectorize(calct_b)

