// Main model in Stan
// Adapted by Samuel Wong from original code by authors of https://github.com/carolinecolijn/distancing-impact-covid19

functions{
    real[] seir(real t,        // time (actual time; not an increment starting at 1)
                real[] state,  // state
                real[] theta,  // parameters
                real[] x_r,    // data (real)
                int[]  x_i) {  // data (integer)
        real S     = state[1];
        real E1    = state[2];
        real E2    = state[3];
        real I     = state[4]; 
        real Q     = state[5];
        real R     = state[6];
        real Sd    = state[7]; 
        real E1d   = state[8];
        real E2d   = state[9];
        real Id    = state[10];
        real Qd    = state[11];
        real Rd    = state[12];
        
        real R0  = theta[1];
        real f2  = theta[2];
        real f3  = theta[3];
        real f4  = theta[4];
        real f5  = theta[5];
        real f6  = theta[6];
        real f7  = theta[7];
        
        real N      = x_r[1];
        real D      = x_r[2];
        real k1     = x_r[3];
        real k2     = x_r[4];
        real q      = x_r[5];
        real r      = x_r[6];
        real ur     = x_r[7];
        real f1     = x_r[8]; 
        real start_decline = x_r[9];
        real end_decline   = x_r[10];
        real fixed_f_forecast = x_r[11];
        real last_day_obs = x_r[12];
        real day_start_fixed_f_forecast = x_r[13];
        
        
        real dydt[12];
        
        real f;
        
        if (t < 14) {
            f = f1;
        }
        if (t >= 14 && t < 21) {
            f = f2 + (21 - t) * (f1 - f2) / 7;
        }
        if (t >= 21 && t < 79) {
            f = f2;
        }
        if (t >= 79 && t < 86) {
            f = f3 + (86 - t) * (f2 - f3) / 7;
        }
        if (t >= 86 && t < 115) {
            f = f3;
        }
        if (t >= 115 && t < 122) {
            f = f4 + (122 - t) * (f3 - f4) / 7;
        }
        if (t >= 122 && t < 196) {
            f = f4;
        }
        if (t >= 196 && t < 203) {
            f = f5 + (203 - t) * (f4 - f5) / 7;
        }
        if (t >= 203 && t < 226) {
            f = f5;
        }
        if (t >= 226 && t < 233) {
            f = f6 + (233 - t) * (f5 - f6) / 7;
        }
        if (t >= 233 && t < 251) {
            f = f6;
        }
        if (t >= 251 && t < 258) {
            f = f7 + (258 - t) * (f6 - f7) / 7;
        }
        if (t >= 258) {
            f = f7;
        }
        if (t >= day_start_fixed_f_forecast && fixed_f_forecast != 0) {
            f = fixed_f_forecast;
        }
        
        
        dydt[1]  = -(R0 /(D +1/k2 )) * (I  + E2  + f*(Id +E2d )) * S /N  - r *S  + ur *Sd ;
        dydt[2]  = (R0 /(D +1/k2 )) * (I  + E2  + f*(Id +E2d )) * S /N  - k1 *E1  -r *E1  + ur *E1d ;
        dydt[3]  = k1 *E1  - k2 *E2  - r *E2  + ur *E2d ;
        dydt[4]  = k2 *E2  - q *I  - I /D  - r *I  + ur *Id ;
        dydt[5]  = q *I  - Q /D  - r *Q  + ur *Qd ;
        dydt[6]  = I /D  + Q /D  - r *R  + ur *Rd ;
        
        dydt[7]  = -(f*R0 /(D +1/k2 )) * (I +E2  + f*(Id +E2d )) * Sd /N  + r *S  - ur *Sd ;
        dydt[8]  = (f*R0 /(D +1/k2 )) * (I +E2  + f*(Id +E2d )) * Sd /N  - k1 *E1d  +r *E1  - ur *E1d ;
        dydt[9]  = k1 *E1d  - k2 *E2d  + r *E2  - ur *E2d ;
        dydt[10] = k2 *E2d  - q *Id  - Id /D  + r *I  - ur *Id ;
        dydt[11] = q *Id  - Qd /D  + r *Q  - ur *Qd ;
        dydt[12] = Id /D  + Qd /D  + r *R  - ur *Rd ;
        
        return dydt;
    }
}

data {
    int<lower=0> T;     // number of time steps
    int<lower=0> N;     // number of days
    int<lower=0> R;     // number of regions
    real y0[R,12];        // initial state
    real t0;            // first time step
    real time[T];       // time increments
    int days[N];        // day increments
    int last_day_obs;   // last day of observed data; days after this are projections
    int daily_cases[R,last_day_obs]; // daily new case counts
    real x_r[R,13];       // data for ODEs (real numbers)
    //real sampFrac[N];   // fraction of cases sampled per time step
    real delayScale;    // Weibull parameter for delay in becoming a case count
    real delayShape;    // Weibull parameter for delay in becoming a case count
    int time_day_id[N]; // last time increment associated with each day
    int time_day_id0[N];// first time increment for Weibull integration of case counts
    real R0_prior[2];   // lognormal log mean and SD for R0 prior
    real phi_prior;     // SD of normal prior on 1/sqrt(phi) [NB2(mu, phi)]
    real f2_prior[2];   // beta prior for f2
    real f3_prior[2];   // beta prior for f3
    real f4_prior[2];   // beta prior for f4
    real f5_prior[2];   // beta prior for f5
    real f6_prior[2];   // beta prior for f6
    real f7_prior[2];   // beta prior for f7
    real psi1_prior[2];   // beta prior for psi1    
    real psi2_prior[2];   // beta prior for psi2
    real psi3_prior[2];   // beta prior for psi3
    real psi4_prior[2];   // beta prior for psi4
    
    int day_inc_sampling;   // day to switch to sampFrac2
    //real sampFrac2_prior[2];   // beta prior for sampFrac2
    
    int<lower=0, upper=1> priors_only; // logical: include likelihood or just priors?
        int<lower=0, upper=1> est_phi; // estimate NB phi?
        //int<lower=0, upper=N> n_sampFrac2; // number of sampFrac2
    int<lower=0, upper=2> obs_model; // observation model: 0 = Poisson, 1 = NB2
    real<lower=0> rw_sigma;
    int tests[N];
    real ode_control[3];
    int<lower=1> N_lik; // number of days in the likelihood
    int dat_in_lik[N_lik]; // vector of data to include in the likelihood
}
transformed data {
    int x_i[0]; // empty; needed for ODE function
}
parameters {
    real R0; // Stan ODE solver seems to be more efficient without this bounded at > 0
    real<lower=0, upper=1> f2[R]; // strength of social distancing
    real<lower=0, upper=1> f3[R]; // strength of social distancing
    real<lower=0, upper=1> f4[R]; // strength of social distancing
    real<lower=0, upper=1> f5[R]; // strength of social distancing
    real<lower=0, upper=1> f6[R]; // strength of social distancing
    real<lower=0, upper=1> f7[R]; // strength of social distancing
    //real<lower=0> phi[R,est_phi]; // NB2 (inverse) dispersion; `est_phi` turns on/off
    real<lower=0> phi[R];
    //real<lower=0, upper=1> sampFrac2[n_sampFrac2];
    real<lower=0, upper=1> sampFrac_1[R];
    real<lower=0, upper=1> sampFrac_2[R];
    real<lower=0, upper=1> sampFrac_3[R];
    real<lower=0, upper=1> sampFrac_4[R];
}
transformed parameters {
    real meanDelay = delayScale * tgamma(1 + 1 / delayShape);
    real dx = time[2] - time[1]; // time increment
    real ft[T];
    real lambda_d[R,N];
    real sum_ft_inner;
    real eta[R,N]; // expected value on link scale (log)
    real k2;
    real E2;
    real E2d;
    real theta[R,7];
    real y_hat[R,T,12];
    real<lower=0> alpha[R,N]; // 1st shape parameter for the beta distribution
    real<lower=0> beta[R,N]; // 2nd shape parameter for the beta distribution
    real<lower=0, upper=1> sampFrac[R,N];
    for (r in 1:R) {
        theta[r,1] = R0;
        theta[r,2] = f2[r];
        theta[r,3] = f3[r];
        theta[r,4] = f4[r];
        theta[r,5] = f5[r];
        theta[r,6] = f6[r];
        theta[r,7] = f7[r];
        
        y_hat[r] = integrate_ode_rk45(seir, y0[r], t0, time, theta[r], x_r[r], x_i, ode_control[1], ode_control[2], ode_control[3]);
        // print("y0: ", y0[r]);
        // print("theta: ", theta[r]);
        // print("y_hat at T: ", y_hat[r,T]);
    
        for (n in 1:N) {
            
            if (n >= 1 && n < 16) {
                sampFrac[r,n] = sampFrac_1[r];
            }
            if (n >= 16 && n < 40) {
                sampFrac[r,n] = sampFrac_2[r];
            }
            if (n >= 40 && n < 53) {
                sampFrac[r,n] = sampFrac_3[r];
            }
            if (n >= 53 && n <= N) {
                sampFrac[r,n] = sampFrac_4[r];
            }

            //print(sampFrac[n]);
            for (t in 1:T) {
                ft[t] = 0; // initialize at 0 across the full 1:T
            }
        
            for (t in time_day_id0[n]:time_day_id[n]) { // t is an increment here
                k2 = x_r[r,4];
                E2 = y_hat[r,t,3];
                E2d = y_hat[r,t,9];
            
                ft[t] = sampFrac[r,n] * k2 * (E2 + E2d) *
                    exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delayShape, delayScale));
            }
            
            sum_ft_inner = 0;
            for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
                sum_ft_inner += ft[t];
            }
            lambda_d[r,n] = 0.5 * dx *
                (ft[time_day_id0[n]] + 2 * sum_ft_inner + ft[time_day_id[n]]);
            eta[r,n] = log(lambda_d[r,n]);
        }
        //print("component ", r, " lambdas: ", lambda_d);
    
        if (obs_model == 2) { // Beta-Binomial
            for (n in 1:N) {
                eta[r,n] = inv_logit(exp(eta[r,n]));
                alpha[r,n] = eta[r,n] * phi[r];
                beta[r,n] = (1 - eta[r,n]) * phi[r];
            }
        } else {
            for (n in 1:N) {
                alpha[r,n] = 0;
                beta[r,n] = 0;
            }
        }        
    }
    

}
model {
    // priors:
        if (est_phi && obs_model == 1) { // NB2
            // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
            // D(expression(1/sqrt(x)), "x"); log(0.5 * x^-0.5/sqrt(x)^2
            for (r in 1:R) {
                                                  1/sqrt(phi[r]) ~ normal(0, phi_prior);
                                                  target += log(0.5) - 1.5 * log(phi[r]); // Jacobian adjustment
            }
        }
    if (est_phi && obs_model == 2) { // Beta-Binomial
        for (r in 1:R) {
            phi[r] ~ normal(0, phi_prior);
        }
    }
    R0 ~ lognormal(R0_prior[1], R0_prior[2]);
    f2 ~ beta(f2_prior[1], f2_prior[2]);
    f3 ~ beta(f3_prior[1], f3_prior[2]);
    f4 ~ beta(f4_prior[1], f4_prior[2]);
    f5 ~ beta(f5_prior[1], f5_prior[2]);
    f6 ~ beta(f6_prior[1], f6_prior[2]);
    f7 ~ beta(f7_prior[1], f7_prior[2]);
    
    sampFrac_1 ~ beta(psi1_prior[1], psi1_prior[2]);
    sampFrac_2 ~ beta(psi2_prior[1], psi2_prior[2]);
    sampFrac_3 ~ beta(psi3_prior[1], psi3_prior[2]);
    sampFrac_4 ~ beta(psi4_prior[1], psi4_prior[2]);
    
    // print("phi: ", phi);
    // print("R0: ", R0);
    // print("f2: ", f2);
    // print("sampFrac_1: ", sampFrac_1);
    
    // data likelihood:
    for (r in 1:R) {
        if (!priors_only) {
            if (obs_model == 0) {
                daily_cases[r,dat_in_lik] ~ poisson_log(eta[r,dat_in_lik]);
            } else if (obs_model == 1) {
                daily_cases[r,dat_in_lik] ~ neg_binomial_2_log(eta[r,dat_in_lik], phi[r]);
                
            //    print(eta[r,dat_in_lik]);
            // } else {
            //     daily_cases[r,dat_in_lik] ~ beta_binomial(tests[dat_in_lik], alpha[r,dat_in_lik], beta[r,dat_in_lik]);
            }
        }
    }
}

generated quantities{
    int y_rep[R,N]; // posterior predictive replicates
    
    for (r in 1:R) {
    for (n in 1:N) {
        if (obs_model == 0) {
            y_rep[r,n] = poisson_log_rng(eta[r,n]);
        } else if (obs_model == 1) {
            y_rep[r,n] = neg_binomial_2_log_rng(eta[r,n], phi[r]);
        // } else {
        //     y_rep[r,n] = beta_binomial_rng(tests[n], alpha[r,n], beta[r,n]);
        }
    }
    }
}