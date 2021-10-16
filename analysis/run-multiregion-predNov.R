# Run the hierarchical multiregion model on the real data
# Fit data up to Nov 7, predict to Nov 30

seeiqr_model <- rstan::stan_model("analysis/seeiqr-multiregion-predNov.stan")
source("analysis/fit_seeiqr-multiregion-predNov.R")

# Case counts
Island <- read.csv("data-generated/daily-cases-Island.csv")
Interior <- read.csv("data-generated/daily-cases-Interior.csv")
Coastal <- read.csv("data-generated/daily-cases-Coastal.csv")
Northern <- read.csv("data-generated/daily-cases-Northern.csv")
Fraser <- read.csv("data-generated/daily-cases-Fraser.csv")

daily_diffs <- rbind(Island$cases,Coastal$cases,Northern$cases,Interior$cases,Fraser$cases)[,1:251]
pop_size <- c(843375,1225195,297570,795116,1889225) # regional populations

# priors, given in MODE/sd
f2_prior = c(0.4, 0.25)
f3_prior = c(0.5, 0.25)
f4_prior = c(0.6, 0.25)
f5_prior = c(0.7, 0.25)
f6_prior = c(0.8, 0.25)
#f7_prior = c(0.6, 0.25)  No f7 for this period
psi1pi = c(0.1, 0.2)
psi2pi = c(0.2, 0.2)
psi3pi = c(0.3, 0.2)
psi4pi = c(0.4, 0.2) 

options(mc.cores = parallel::detectCores())
m7 <- fit_seeiqr(daily_diffs,
                 pop_size = pop_size,
                 seeiqr_model = seeiqr_model, iter = 2000, refresh = 1, warmup = 1000, chains = 4,
                 f2_prior = f2_prior,
                 f3_prior = f3_prior,
                 f4_prior = f4_prior,
                 f5_prior = f5_prior,
                 f6_prior = f6_prior,
                 #f7_prior = f7_prior,
                 psi1pi = psi1pi,
                 psi2pi = psi2pi,
                 psi3pi = psi3pi,
                 psi4pi = psi4pi,
                 save_state_predictions = TRUE
)


save(f2_prior, f3_prior, f4_prior, f5_prior, f6_prior, psi1pi, psi2pi, psi3pi, psi4pi, m7, file="data-generated/fit-predNov-multiregion.rda")

