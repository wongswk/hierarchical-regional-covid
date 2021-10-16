#' Fit the Stan SEEIQR hierarchical multiregion model
#' Adapted by Samuel Wong from original code by authors of https://github.com/carolinecolijn/distancing-impact-covid19
#'
#' @param daily_cases Vectors of daily new cases
#' @param daily_tests An optional vector of daily test numbers. Should include
#'   assumed tests for the forecast. I.e. `nrow(daily_cases) + forecast_days =
#'   length(daily_tests)`. Only used in the case of the beta-binomial (which
#'   isn't working very well).
#' @param Model from `rstan::stan_model(seeiqr_model)`.
#' @param obs_model Type of observation model
#' @param forecast_days Number of days into the future to forecast. The model
#'   will run slightly faster with fewer forecasted days.
#' @param time_increment Time increment for ODEs and Weibull delay-model
#'   integration
#' @param days_back Number of days to go back for Weibull delay-model
#'   integration
#' @param R0_prior Lognormal log mean and SD for R0 prior
#' @param phi_prior SD of `1/sqrt(phi) ~ Normal(0, SD)` prior, where NB2(mu,
#'   phi) and `Var(Y) = mu + mu^2 / phi`.
#'   <https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations>
#' @param f2_prior Beta mean and SD for `f2` parameter
#' @param f3_prior Beta mean and SD for `f3` parameter
#' @param f4_prior Beta mean and SD for `f4` parameter
#' @param f5_prior Beta mean and SD for `f5` parameter
#' @param f6_prior Beta mean and SD for `f6` parameter
#' @param f7_prior Beta mean and SD for `f7` parameter 
#' @param sampFrac2_prior `sampFrac` prior starting on
#'   `sampled_fraction_day_change` if `sampFrac2_type` is "estimated" or "rw".
#'   In the case of the random walk, this specifies the initial state prior. The
#'   two values correspond to the mean and SD of a Beta distribution.
#' @param sampFrac2_type How to treat the sample fraction. Fixed, estimated, or
#'   a constrained random walk.
#' @param rw_sigma The standard deviation on the sampFrac2 random walk.
#' @param seed MCMC seed
#' @param chains Number of MCMC chains
#' @param iter MCMC iterations per chain
#' @param sampled_fraction1 Fraction sampled before
#'   `sampled_fraction_day_change`
#' @param sampled_fraction2 Fraction sampled after
#'   `sampled_fraction_day_change` 
#' @param sampled_fraction_day_change Date fraction sample changes
#' @param sampled_fraction_vec An optional vector of sampled fractions. Should
#'   be of length: `nrow(daily_cases) + forecast_days`. 
#' @param fixed_f_forecast Optional fixed `f` for forecast.
#' @param pars A named numeric vector of fixed parameter values
#' @param i0 A scaling factor FIXME
#' @param fsi Fraction socially distancing. Derived parameter.
#' @param nsi Fraction not socially distancing. Derived parameter.
#' @param state_0 Initial state: a named numeric vector
#' @param save_state_predictions Include the state predictions? `y_hat` Will
#'   make the resulting model object much larger.
#' @param delayScale Weibull scale parameter for the delay in reporting.
#' @param delayShape Weibull shape parameter for the delay in reporting.
#' @param ode_control Control options for the Stan ODE solver. First is relative
#'   difference, that absolute difference, and then maximum iterations. The
#'   values here are the Stan defaults.
#' @param daily_cases_omit An optional vector of days to omit from the data
#'   likelihood.
#' @param ... Other arguments to pass to [rstan::sampling()].
#' @author Sean Anderson

fit_seeiqr <- function(daily_cases,  # will be a matrix
                       daily_tests = NULL,
                       seeiqr_model,
                       obs_model = c("NB2", "Poisson", "beta-binomial"),
                       forecast_days = 1, 
                       time_increment = 0.1,
                       days_back = 45,
                       R0_prior = c(log(2.6), 0.2),
                       phi_prior = 1,
                       f2_prior,
                       f3_prior,
                       f4_prior,
                       f5_prior,
                       f6_prior,
                       f7_prior,
                       psi1pi,
                       psi2pi,
                       psi3pi,
                       psi4pi,
                       #sampFrac2_prior = c(0.4, 0.2),
                       #sampFrac2_type = c("fixed", "estimated", "rw"),
                       rw_sigma = 0.1,
                       seed = 42,
                       chains = 8,
                       iter = 1000,
                       #sampled_fraction1 = 0.1,
                       #sampled_fraction2 = 0.3,
                       sampled_fraction_day_change = 14,
                       sampled_fraction_vec = NULL,
                       fixed_f_forecast = NULL,
                       day_start_fixed_f_forecast = ncol(daily_cases) + 1,
                       pars = c(
                         N = 5.1e6, D = 5, k1 = 1 / 5, # note, we fill N with the actual region population later
                         k2 = 1, q = 0.05,
                         r = 0.1, ur = 0.02, f1 = 1.0,
                         start_decline = 15,
                         end_decline = 22
                       ),
                       i0 = 8,
                       fsi = pars[["r"]] / (pars[["r"]] + pars[["ur"]]),
                       nsi = 1 - fsi,
                       # state_0 = c(
                       #   S = nsi * (pars[["N"]] - i0),
                       #   E1 = 0.4 * nsi * i0,
                       #   E2 = 0.1 * nsi * i0,
                       #   I = 0.5 * nsi * i0,
                       #   Q = 0,
                       #   R = 0,
                       #   Sd = fsi * (pars[["N"]] - i0),
                       #   E1d = 0.4 * fsi * i0,
                       #   E2d = 0.1 * fsi * i0,
                       #   Id = 0.5 * fsi * i0,
                       #   Qd = 0,
                       #   Rd = 0
                       # ),
                       pop_size,
                       save_state_predictions = FALSE,
                       delayScale = 9.85,
                       delayShape = 1.73,
                       ode_control = c(1e-6, 1e-6, 1e6),
                       daily_cases_omit = NULL,
                       ...) {
  obs_model <- match.arg(obs_model)
  obs_model <-
    if (obs_model == "Poisson") {
      0L
    } else if (obs_model == "NB2") {
      1L
    } else {
      2L
    }

  R <- nrow(daily_cases)
  #x_r <- matrix(rep(pars,R), nrow=R)
  x_r <- pars
  pop_total <- sum(pop_size)
  
  state_0 <- matrix(NA, nrow=R, ncol = 12)
  for (r in 1:R) {
    state_0[r,] = c(
      S = nsi * (pop_size[r] - i0 * pop_size[r]/pop_total),
      E1 = 0.4 * nsi * i0* pop_size[r]/pop_total,
      E2 = 0.1 * nsi * i0* pop_size[r]/pop_total,
      I = 0.5 * nsi * i0* pop_size[r]/pop_total,
      Q = 0,
      R = 0,
      Sd = fsi * (pop_size[r] - i0* pop_size[r]/pop_total),
      E1d = 0.4 * fsi * i0* pop_size[r]/pop_total,
      E2d = 0.1 * fsi * i0* pop_size[r]/pop_total,
      Id = 0.5 * fsi * i0* pop_size[r]/pop_total,
      Qd = 0,
      Rd = 0
    )
    #x_r[r,1] = pop_size[r]
  }
    
  
  
  # if (!is.null(daily_tests)) {
  #   stopifnot(length(daily_cases) + forecast_days == length(daily_tests))
  #   if (min(daily_tests) == 0) {
  #     warning("Replacing 0 daily tests with 1.")
  #     daily_tests[daily_tests == 0] <- 1
  #   }
  # }
  
  stopifnot(
    names(x_r) ==
      c("N", "D", "k1", "k2", "q", "r", "ur", "f1", "start_decline", "end_decline")
  )
  stopifnot(
    names(state_0) == c("S", "E1", "E2", "I", "Q", "R", "Sd", "E1d", "E2d", "Id", "Qd", "Rd")
  )
  
  days <- seq(1, ncol(daily_cases) + forecast_days)
  last_day_obs <- ncol(daily_cases)
  time <- seq(-30, max(days), time_increment)
  x_r <- c(x_r, if (!is.null(fixed_f_forecast)) fixed_f_forecast else 0)
  names(x_r)[length(x_r)] <- "fixed_f_forecast"
  x_r <- c(x_r, c("last_day_obs" = last_day_obs))
  x_r <- c(x_r, c("day_start_fixed_f_forecast" = day_start_fixed_f_forecast))
  
  x_r <- matrix(rep(x_r,each=R), nrow=R)
  for (r in 1:R)
    x_r[r,1] = pop_size[r]
  
  # find the equivalent time of each day (end):
  get_time_id <- function(day, time) max(which(time <= day))
  time_day_id <- vapply(days, get_time_id, numeric(1), time = time)
  
  get_time_day_id0 <- function(day, time, days_back) {
    # go back `days_back` or to beginning if that's negative time:
    check <- time <= (day - days_back)
    if (sum(check) == 0L) {
      1L
    } else {
      max(which(check))
    }
  }
  # find the equivalent time of each day (start):
  time_day_id0 <- vapply(days, get_time_day_id0, numeric(1),
                         time = time, days_back = days_back
  )

  
  beta_sd <- f2_prior[2]
  beta_mean <- f2_prior[1]
  beta_sdf3 <- f3_prior[2]
  beta_meanf3 <- f3_prior[1]
  beta_sdf4 <- f4_prior[2]
  beta_meanf4 <- f4_prior[1]
  beta_sdf5 <- f5_prior[2]
  beta_meanf5 <- f5_prior[1]
  beta_sdf6 <- f6_prior[2]
  beta_meanf6 <- f6_prior[1]
  beta_sdf7 <- f7_prior[2]
  beta_meanf7 <- f7_prior[1] 
  beta_shape1 <- get_beta_params(beta_mean, beta_sd)$alpha
  beta_shape2 <- get_beta_params(beta_mean, beta_sd)$beta
  beta_shapef31 <- get_beta_params(beta_meanf3, beta_sdf3)$alpha
  beta_shapef32 <- get_beta_params(beta_meanf3, beta_sdf3)$beta
  beta_shapef41 <- get_beta_params(beta_meanf4, beta_sdf4)$alpha
  beta_shapef42 <- get_beta_params(beta_meanf4, beta_sdf4)$beta
  beta_shapef51 <- get_beta_params(beta_meanf5, beta_sdf5)$alpha
  beta_shapef52 <- get_beta_params(beta_meanf5, beta_sdf5)$beta
  beta_shapef61 <- get_beta_params(beta_meanf6, beta_sdf6)$alpha
  beta_shapef62 <- get_beta_params(beta_meanf6, beta_sdf6)$beta
  beta_shapef71 <- get_beta_params(beta_meanf7, beta_sdf7)$alpha
  beta_shapef72 <- get_beta_params(beta_meanf7, beta_sdf7)$beta
  
  bpsi1 <- get_beta_params(psi1pi[1], psi1pi[2])
  bpsi2 <- get_beta_params(psi2pi[1], psi2pi[2])
  bpsi3 <- get_beta_params(psi3pi[1], psi3pi[2])
  bpsi4 <- get_beta_params(psi4pi[1], psi4pi[2])
  
  dat_in_lik <- seq(1, last_day_obs)
  if (!is.null(daily_cases_omit)) {
    dat_in_lik <- dat_in_lik[-daily_cases_omit]
  }
  
  #show(daily_cases)
  show(state_0)
  show(x_r)
  
  stan_data <- list(
    R = R,
    T = length(time),
    days = days,
    daily_cases = daily_cases,
    tests = if (is.null(daily_tests)) rep(log(1), length(days)) else daily_tests,
    N = length(days),
    y0 = state_0,
    t0 = min(time) - 0.000001,
    time = time,
    x_r = x_r,
    delayShape = delayShape,
    delayScale = delayScale,
    #sampFrac = sampFrac,
    time_day_id = time_day_id,
    time_day_id0 = time_day_id0,
    R0_prior = R0_prior,
    phi_prior = phi_prior,
    f2_prior = c(beta_shape1, beta_shape2),
    f3_prior = c(beta_shapef31, beta_shapef32),
    f4_prior = c(beta_shapef41, beta_shapef42),
    f5_prior = c(beta_shapef51, beta_shapef52),
    f6_prior = c(beta_shapef61, beta_shapef62),
    f7_prior = c(beta_shapef71, beta_shapef72),
    psi1_prior = c(bpsi1$alpha, bpsi1$beta),
    psi2_prior = c(bpsi2$alpha, bpsi2$beta),
    psi3_prior = c(bpsi3$alpha, bpsi3$beta),
    psi4_prior = c(bpsi4$alpha, bpsi4$beta),
    #sampFrac2_prior = c(sampFrac2_beta_shape1, sampFrac2_beta_shape2),
    day_inc_sampling = sampled_fraction_day_change,
    #n_sampFrac2 = n_sampFrac2,
    rw_sigma = rw_sigma,
    priors_only = 0L,
    last_day_obs = last_day_obs,
    obs_model = obs_model,
    ode_control = ode_control,
    est_phi = if (obs_model %in% c(1L, 2L)) 1L else 0L,
    dat_in_lik = dat_in_lik,
    N_lik = length(dat_in_lik)
  )

  initf <- function(stan_data) {
    # R0 <- rlnorm(1, R0_prior[1], R0_prior[2])
    # f2 <- rbeta(
    #   nrow(daily_cases),
    #   get_beta_params(f2_prior[1], f2_prior[2])$alpha,
    #   get_beta_params(f2_prior[1], f2_prior[2])$beta
    # )
    # f3 <- rbeta(
    #   nrow(daily_cases),
    #   get_beta_params(f3_prior[1], f3_prior[2])$alpha,
    #   get_beta_params(f3_prior[1], f3_prior[2])$beta    
    # )
    # f4 <- rbeta(
    #   nrow(daily_cases),
    #   get_beta_params(f4_prior[1], f4_prior[2])$alpha,
    #   get_beta_params(f4_prior[1], f4_prior[2])$beta    
    # )
    # f5 <- rbeta(
    #   nrow(daily_cases),
    #   get_beta_params(f5_prior[1], f5_prior[2])$alpha,
    #   get_beta_params(f5_prior[1], f5_prior[2])$beta    
    # )
    # f6 <- rbeta(
    #   nrow(daily_cases),
    #   get_beta_params(f6_prior[1], f6_prior[2])$alpha,
    #   get_beta_params(f6_prior[1], f6_prior[2])$beta    
    # )
    # 
    # sampFrac_1 <- rbeta(nrow(daily_cases), 3.5, 31.5) #0.1, sd=0.05
    # sampFrac_2 <- rbeta(nrow(daily_cases), 24.9, 58.1) #0.3, sd=0.05
    # sampFrac_3 <- rbeta(nrow(daily_cases), 38, 57) #0.4, sd=0.05
    # sampFrac_4 <- rbeta(nrow(daily_cases), 49.5, 49.5) #0.5, sd=0.05 
    
    R0 <- 2.6
    f2 <- c(0.45,0.65,0.35,0.6,0.5)
    f3 <- c(0.6,0.3,0.75,0.4,0.6)
    f4 <- c(0.8,0.75,0.75,0.9,0.69)
    f5 <- c(0.75,0.8,0.6,0.8,0.67)
    f6 <- c(0.8,0.7,0.64,0.85,0.71)
    f7 <- rep(0.6,5)
    phi <- c(3,2,8,3,8)
    # sampFrac_1 <- rep(0.1,R)
    # sampFrac_2 <- rep(0.2,R)
    # sampFrac_3 <- rep(0.3,R)
    # sampFrac_4 <- rep(0.4,R)
    sampFrac_1 <- c(0.1, 0.3, 0.05, 0.05, 0.3)
    sampFrac_2 <- c(0.18, 0.39, 0.15, 0.24, 0.36)
    sampFrac_3 <- c(0.41, 0.46, 0.37, 0.4, 0.47)
    sampFrac_4 <- rep(0.5, 5)
    
    f2 <- as.array(f2[1:R])
    f3 <- as.array(f3[1:R])
    f4 <- as.array(f4[1:R])
    f5 <- as.array(f5[1:R])
    f6 <- as.array(f6[1:R])
    f7 <- as.array(f7[1:R])
    phi <- as.array(phi[1:R])
    sampFrac_1 <- as.array(sampFrac_1[1:R])
    sampFrac_2 <- as.array(sampFrac_2[1:R])
    sampFrac_3 <- as.array(sampFrac_3[1:R])
    sampFrac_4 <- as.array(sampFrac_4[1:R])
    
    
    init <- list(R0 = R0, f2 = f2, f3 = f3, f4 = f4, f5 = f5, f6 = f6, f7 = f7, phi = phi, sampFrac_1 = sampFrac_1, sampFrac_2 = sampFrac_2, sampFrac_3 = sampFrac_3, sampFrac_4 = sampFrac_4)
    init
  }
  pars_save <- c("R0", "f2", "f3", "f4", "f5", "f6", "f7", "phi", "lambda_d",  "sampFrac_1", "sampFrac_2", "sampFrac_3", "sampFrac_4", "y_rep")###
  if (save_state_predictions) pars_save <- c(pars_save, "y_hat", "sampFrac") ###
  fit <- rstan::sampling(
    seeiqr_model,
    data = stan_data,
    iter = iter,
    chains = chains,
    init = function() initf(stan_data),
    seed = seed, # https://xkcd.com/221/
    pars = pars_save,
    ... = ...
  )
  post <- rstan::extract(fit)
  list(
    fit = fit, post = post, phi_prior = phi_prior, R0_prior = R0_prior,
    f2_prior = f2_prior, f3_prior = f3_prior, f4_prior = f4_prior, f5_prior = f5_prior, f6_prior = f6_prior, f7_prior = f7_prior, obs_model = obs_model, state_0 = state_0,
    daily_cases = daily_cases, daily_tests = daily_tests, days = days, time = time,
    last_day_obs = last_day_obs, pars = x_r, f2_prior_beta_shape1 = beta_shape1,
    f2_prior_beta_shape2 = beta_shape2, f3_prior_beta_shape1 = beta_shapef31,
    f3_prior_beta_shape2 = beta_shapef32, f4_prior_beta_shape1 = beta_shapef41,
    f4_prior_beta_shape2 = beta_shapef42, f5_prior_beta_shape1 = beta_shapef51,
    f5_prior_beta_shape2 = beta_shapef52, f6_prior_beta_shape1 = beta_shapef61,
    f6_prior_beta_shape2 = beta_shapef62, f7_prior_beta_shape1 = beta_shapef71, f7_prior_beta_shape2 = beta_shapef72, stan_data = stan_data
  )
}

# get_beta_params <- function(mu, sd) {
#   var <- sd^2
#   alpha <- ((1 - mu) / var - 1 / mu) * mu^2
#   beta <- alpha * (1 / mu - 1)
#   list(alpha = alpha, beta = beta)
# }

get_beta_params <- function(mode, sd) {
  var <- sd^2
  
  f <- function(alpha, mode, sd) {
    beta <- (alpha-1)/mode + 2 - alpha
    
    sd^2 - (alpha * beta)/ ((alpha+beta)^2 * (alpha + beta+1))
  }
  
  foo <- uniroot(f, c(1,100), mode, sd)
  
  alpha <- foo$root
  beta <- (alpha-1)/mode + 2 - alpha
  list(alpha = alpha, beta = beta)
}
