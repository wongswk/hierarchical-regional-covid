# Generate simulated case counts using parameters in simulation-params.R
# Adapted by Samuel Wong from original code by authors of https://github.com/carolinecolijn/distancing-impact-covid19

#' Lambda_d: returns expected number of cases on day d
#' @author Caroline Colijn, Jessica Stockdale
getlambd <- function(out,
                     pars,
                     day,
                     #data = bcdata,
                     #sampFrac = 0.1,
                     delayShape = 1.73,
                     delayScale = 9.85) {
  meanDelay <- delayScale * gamma(1 + 1 / delayShape)
  try(if (var(diff(out$time)) > 0.005) {
    stop("approx integral assumes equal time steps")
  })
  try(if (max(out$time) < day) {
    stop("model simulation is not long enough for the data")
  })
  try(if (min(out$time) > day - (2 * meanDelay + 1)) {
    stop("we need an earlier start time for the model")
  })
  # relevant times to identify new cases
  ii <- which(out$time > day - 45 & out$time <= day)
  dx <- out$time[ii[2]] - out$time[ii[1]] # assumes equal time intervals
  # all new cases arising at each of those times
  incoming <- with(pars, {
    k2 * (out$E2[ii] + out$E2d[ii])
  })

  # sampfrac is needed here
  if (day >= 1 & day < 16) {
    thisSamp = pars$sampFrac_1;
  }
  if (day >= 16 & day < 40) {
    thisSamp = pars$sampFrac_2;
  }
  if (day >= 40 & day < 53) {
    thisSamp = pars$sampFrac_3;
  }
  if (day >= 53) {
    thisSamp = pars$sampFrac_4;
  }  
  
  
  # each of the past times' contribution to this day's case count
  ft <- thisSamp * incoming * dweibull(
    x = max(out$time[ii]) - out$time[ii],
    shape = delayShape,
    scale = delayScale
  )
  # return numerical integral of ft
  return(0.5 * (dx) * (ft[1] + 2 * sum(ft[2:(length(ft) - 1)]) + ft[length(ft)]))
}

#' @title Social Distancing Model
#' @author Caroline Colijn
#' @description SEIR-type model with time-dependent social distancing. Social
#'   distancing reduces frequency of contact. Individuals can move between
#'   distanced and not distanced compartments.
#' @param t time
#' @param state (S, E1, E2, I, Q, R, Sd, E1d, E2d, Id, Qd, Rd) S: Susceptible,
#'   E1: Exposed but not infectious, E2: Exposed and Infectious, I: Infectious,
#'   can be quarantined, R: Removed. The d compartments denote socially
#'   distanced individuals.
#' @param pars (N, D, R0, k1, k2, q, r, ur, f) f: strength of social distancing,
#'   r/(r+ur): frac of population who are distancing
#' @param sdtiming timing of social distancing
#' @return time derivatives for input to ODE solver
socdistmodel <- function(t,state,parms) {
  with(as.list(c(
    state, parms
  )), {
    #f <- sdtiming(t, f1 = parms$f1, f2 = parms$f2)
    if (t < 14) {
      f = parms$f1;
    }
    if (t >= 14 & t < 21) {
      f = parms$f2 + (21 - t) * (parms$f1 - parms$f2) / 7;
    }
    if (t >= 21 & t < 79) {
      f = parms$f2;
    }
    if (t >= 79 & t < 86) {
      f = parms$f3 + (86 - t) * (parms$f2 - parms$f3) / 7;
    }
    if (t >= 86 & t < 115) {
      f = parms$f3;
    }
    if (t >= 115 & t < 122) {
      f = parms$f4 + (122 - t) * (parms$f3 - parms$f4) / 7;
    }
    if (t >= 122 & t < 196) {
      f = parms$f4;
    }
    if (t >= 196 & t < 203) {
      f = parms$f5 + (203 - t) * (parms$f4 - parms$f5) / 7;
    }
    if (t >= 203 & t < 226) {
      f = parms$f5;
    }
    if (t >= 226 & t < 233) {
      f = parms$f6 + (233 - t) * (parms$f5 - parms$f6) / 7;
    }
    if (t >= 233 & t < 251) {
      f = parms$f6;
    }
    if (t >= 251 & t < 258) {
      f = parms$f7 + (258 - t) * (parms$f6 - parms$f7) / 7;
    }
    if (t >= 258) {
      f = parms$f7;
    }
    
    dSdt <- -(R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) * S / N - r * S + ur * Sd
    dE1dt <- (R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) * S / N - k1 * E1 - r * E1 + ur * E1d
    dE2dt <- k1 * E1 - k2 * E2 - r * E2 + ur * E2d
    dIdt <- k2 * E2 - q * I - I / D - r * I + ur * Id
    dQdt <- q * I - Q / D - r * Q + ur * Qd
    dRdt <- I / D + Q / D - r * R + ur * Rd
    
    dSddt <- -(f * R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) * Sd / N + r * S - ur * Sd
    dE1ddt <- (f * R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) * Sd / N - k1 * E1d + r * E1 - ur * E1d
    dE2ddt <- k1 * E1d - k2 * E2d + r * E2 - ur * E2d
    dIddt <- k2 * E2d - q * Id - Id / D + r * I - ur * Id
    dQddt <- q * Id - Qd / D + r * Q - ur * Qd
    dRddt <- Id / D + Qd / D + r * R - ur * Rd
    list(c(
      dSdt,
      dE1dt,
      dE2dt,
      dIdt,
      dQdt,
      dRdt,
      dSddt,
      dE1ddt,
      dE2ddt,
      dIddt,
      dQddt,
      dRddt
    ))
  })
}

R <- 5
# load param values here
source("analysis/simulation-params.R")
pop_total <- sum(pop_size)

sim_dat <- list()
ode_out <- list()


for (r in 1:R) {
  pars <- list(
    N = pop_size[r],
    D = 5,
    #other params are as in fit_seeiqr 93-97
    R0 = R0b,
    # ours
    k1 = 1 / 5,
    k2 = 1,
    q = 0.05,
    r = 0.1,
    ur = 0.02,
    sampFrac_1 = sampFrac_1[r],
    sampFrac_2 = sampFrac_2[r],
    sampFrac_3 = sampFrac_3[r],
    sampFrac_4 = sampFrac_4[r],
    f1 = f1[r],
    f2 = f2[r],
    f3 = f3[r],
    f4 = f4[r],
    f5 = f5[r],
    f6 = f6[r],
    f7 = f7[r]
  )
  fsi <- with(pars,
              r / (r + ur))
  nsi <- 1 - fsi
  i0 <- 8
  # state_0 <- c(
  #   S = nsi * (pars$N - i0),
  #   E1 = 0.4 * nsi * i0,
  #   E2 = 0.1 * nsi * i0,
  #   I = 0.5 * nsi * i0,
  #   Q = 0,
  #   R = 0,
  #   Sd = fsi * (pars$N - i0),
  #   E1d = 0.4 * fsi * i0,
  #   E2d = 0.1 * fsi * i0,
  #   Id = 0.5 * fsi * i0,
  #   Qd = 0,
  #   Rd = 0
  state_0 <- c(
    S = nsi * (pop_size[r] - i0 * pop_size[r] / pop_total),
    E1 = 0.4 * nsi * i0 * pop_size[r] / pop_total,
    E2 = 0.1 * nsi * i0 * pop_size[r] / pop_total,
    I = 0.5 * nsi * i0 * pop_size[r] / pop_total,
    Q = 0,
    R = 0,
    Sd = fsi * (pop_size[r] - i0 * pop_size[r] / pop_total),
    E1d = 0.4 * fsi * i0 * pop_size[r] / pop_total,
    E2d = 0.1 * fsi * i0 * pop_size[r] / pop_total,
    Id = 0.5 * fsi * i0 * pop_size[r] / pop_total,
    Qd = 0,
    Rd = 0
  )
  
  times <- seq(from = -30,
               to = 306,
               by = 0.1)
  
  set.seed(seed)
  
  example_simulation <- as.data.frame(deSolve::ode(
    y = state_0,
    times = times,
    func = socdistmodel,
    parms = pars,
  ))
  ode_out[[r]] <- example_simulation
  
  dat <- data.frame(Date = seq(
    lubridate::ymd("2020-03-01"),
    lubridate::ymd("2020-03-01") + max(times) - 1,
    by = "day"
  ))
  dat$day <- seq_along(dat$Date)
  lambda_d <- sapply(seq(1, max(times)), function(x) {
    getlambd(example_simulation, pars = pars, day = x)
  })
  
  sim_datr <- data.frame(
    day = seq(1, max(times)),
    lambda_d = lambda_d,
    obs = MASS::rnegbin(max(times), lambda_d, theta = theta[r])
  )
  sim_dat[[r]] <- sim_datr
  
}

# Check if simulated case counts are reasonable:
# sum(sim_dat[[1]]$obs)
# sum(Island$cases)
# sum(sim_dat[[2]]$obs)
# sum(Coastal$cases)
# sum(sim_dat[[3]]$obs)
# sum(Northern$cases)
# sum(sim_dat[[4]]$obs)
# sum(Interior$cases)
# sum(sim_dat[[5]]$obs)
# sum(Fraser$cases)
