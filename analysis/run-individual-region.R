# Fit models to individual regions on the real data

seeiqr_model <- rstan::stan_model("analysis/seeiqr-multiregion.stan")
source("analysis/fit_seeiqr-multiregion.R")

Island <- read.csv("data-generated/daily-cases-Island.csv")
Interior <- read.csv("data-generated/daily-cases-Interior.csv")
Coastal <- read.csv("data-generated/daily-cases-Coastal.csv")
Northern <- read.csv("data-generated/daily-cases-Northern.csv")
Fraser <- read.csv("data-generated/daily-cases-Fraser.csv")

daily_diffs_each <- rbind(Island$cases,Coastal$cases,Northern$cases,Interior$cases,Fraser$cases)
pop_size_each <- c(843375,1225195,297570,795116,1889225)

mIndiv <- list()

# priors, given in MODE/sd
f2_prior = c(0.4, 0.25)
f3_prior = c(0.5, 0.25)
f4_prior = c(0.6, 0.25)
f5_prior = c(0.7, 0.25)
f6_prior = c(0.8, 0.25)
f7_prior = c(0.6, 0.25)
psi1pi = c(0.1, 0.2)
psi2pi = c(0.2, 0.2)
psi3pi = c(0.3, 0.2)
psi4pi = c(0.4, 0.2) 

for (i in 1:5) {
  show(i)
  pop_size <- pop_size_each[i]
  daily_diffs <- as.matrix(t(daily_diffs_each[i,]))
  show(pop_size)
  show(daily_diffs)
  
  options(mc.cores = parallel::detectCores())
  mIndiv[[i]] <- fit_seeiqr(daily_diffs,
                   pop_size = pop_size,
                   seeiqr_model = seeiqr_model, iter = 2000, refresh = 1, warmup = 1000, chains = 4,
                   f2_prior = f2_prior,
                   f3_prior = f3_prior,
                   f4_prior = f4_prior,
                   f5_prior = f5_prior,
                   f6_prior = f6_prior,
                   f7_prior = f7_prior,
                   psi1pi = psi1pi,
                   psi2pi = psi2pi,
                   psi3pi = psi3pi,
                   psi4pi = psi4pi,                 
                   save_state_predictions = TRUE
  )

}

saveRDS(mIndiv, file="data-generated/fit-endDec-individual-regions.rds")

