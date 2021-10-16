# Make plots based on multiregion and provincial-wide fits
# Needs results saved from "run-multiregion.R" (or "run-multiregion-simulated.R" for simulation study) and "run-provincewide.R"
# Some plotting code adapted by Samuel Wong from original code by authors of https://github.com/carolinecolijn/distancing-impact-covid19

library(rstan)

# Set to TRUE for real data analysis, FALSE for simulation analysis
realdata <- TRUE

# Real data
if (realdata) {
  Island <- read.csv("data-generated/daily-cases-Island.csv")
  Interior <- read.csv("data-generated/daily-cases-Interior.csv")
  Coastal <- read.csv("data-generated/daily-cases-Coastal.csv")
  Northern <- read.csv("data-generated/daily-cases-Northern.csv")
  Fraser <- read.csv("data-generated/daily-cases-Fraser.csv")
  daily_diffs <- rbind(Island$cases,Coastal$cases,Northern$cases,Interior$cases,Fraser$cases)
  
  load("data-generated/fit-endDec-multiregion.rda")
} else {
  load("data-generated/fit-endDec-multiregion-sim-12345.rda")  
}

pop_size <- c(843375,1225195,297570,795116,1889225)
region_names <- c("Island", "Coastal", "Northern", "Interior", "Fraser")

# Provincial results
mprov <- readRDS("data-generated/fit-endDec-provincewide.rds")

library(magrittr)
library(dplyr) 
library(ggplot2)

variables_df <- dplyr::tibble(
  #variable = names(obj$state_0),
  variable = c("S", "E1", "E2", "I", "Q", "R", "Sd", "E1d", "E2d", "Id", "Qd", "Rd"),
  variable_num = 1:12 # seq_along(obj$state_0)
)

.hist_blue <- RColorBrewer::brewer.pal(6, "Blues")[5]
.start <- lubridate::ymd_hms("2020-03-01 00:00:00")
ts_df <- dplyr::tibble(time = m7$time, time_num = seq_along(m7$time))




statesp <- reshape2::melt(mprov$post$y_hat) %>%
  #dplyr::rename(time_num = Var2, variable_num = Var3) %>%
  dplyr::rename(time_num = Var3, variable_num = Var4) %>%
  dplyr::left_join(variables_df, by = "variable_num") %>%
  dplyr::left_join(ts_df, by = "time_num")
prevalencep <- statesp %>%
  dplyr::filter(variable %in% c("I", "Id")) %>%
  group_by(iterations, time) %>%
  summarize(
    I = value[variable == "I"], Id = value[variable == "Id"],
    prevalence = I + Id
  ) %>%
  mutate(day = .start + lubridate::ddays(time))
tempp <- prevalencep %>%
  dplyr::filter(time < ncol(daily_diffs)) %>%
  group_by(day) %>%
  summarise(
    lwr = quantile(prevalence, probs = 0.05),
    lwr2 = quantile(prevalence, probs = 0.25),
    upr = quantile(prevalence, probs = 0.95),
    upr2 = quantile(prevalence, probs = 0.75),
    med = median(prevalence)
  )

for (region in 1:5) {
  
  #  ts_df <- dplyr::tibble(time = m7$time, time_num = seq_along(m7$time))
  states <- reshape2::melt(m7$post$y_hat) %>%
    dplyr::rename(time_num = Var3, variable_num = Var4) %>%
    #dplyr::filter(iterations %in% draws) %>%
    dplyr::filter(Var2 %in% region) %>%
    dplyr::left_join(variables_df, by = "variable_num") %>%
    dplyr::left_join(ts_df, by = "time_num")
  
  prevalence <- states %>%
    dplyr::filter(variable %in% c("I", "Id")) %>%
    group_by(iterations, time) %>%
    summarize(
      I = value[variable == "I"], Id = value[variable == "Id"],
      prevalence = I + Id
    ) %>%
    mutate(day = .start + lubridate::ddays(time))
  
  temp <- prevalence %>%
    dplyr::filter(time < ncol(daily_diffs)) %>%
    group_by(day) %>%
    summarise(
      lwr = quantile(prevalence, probs = 0.05),
      lwr2 = quantile(prevalence, probs = 0.25),
      upr = quantile(prevalence, probs = 0.95),
      upr2 = quantile(prevalence, probs = 0.75),
      med = median(prevalence)
    )
  
  g_prev <- ggplot(temp, aes(x = day, y = med, ymin = lwr, ymax = upr, colour = .hist_blue, fill = .hist_blue)) +
    ylab("Modelled prevalence") +
    ggtitle(region_names[region]) + 
    coord_cartesian(expand = FALSE, xlim = c(.start, .start + lubridate::ddays(m7$last_day_obs + 1)), ylim = c(0, max(c(temp$upr, tempp$upr * pop_size[region] / sum(pop_size))) * 1.04)) +
    xlab("") + 
    geom_line(data = temp, aes(x = day, y = med, colour = 1), alpha = 1, lwd = 1, inherit.aes = FALSE, colour=.hist_blue) +
    geom_line(data = tempp, aes(x = day, y = med * pop_size[region] / sum(pop_size), colour = 1), alpha = 1, lwd = 1, inherit.aes = FALSE, colour="dark red") +
    geom_ribbon(data = tempp, alpha = 0.2, fill = "dark red", colour = NA, mapping = aes(ymin = lwr * pop_size[region] / sum(pop_size), ymax = upr* pop_size[region] / sum(pop_size))) +
    geom_ribbon(data = tempp, alpha = 0.2, mapping = aes(ymin = lwr2* pop_size[region] / sum(pop_size), ymax = upr2* pop_size[region] / sum(pop_size)), fill = "dark red", colour = NA) +
    geom_ribbon(data = temp, alpha = 0.2, fill = .hist_blue, colour = NA, mapping = aes(ymin = lwr, ymax = upr)) +
    geom_ribbon(data = temp, alpha = 0.2, mapping = aes(ymin = lwr2, ymax = upr2), fill = .hist_blue, colour = NA) +
    geom_vline(xintercept = .start + lubridate::ddays(14), lty = 2, col = "grey40") +
    geom_vline(xintercept = .start + lubridate::ddays(79), lty = 2, col = "grey40") +
    geom_vline(xintercept = .start + lubridate::ddays(115), lty = 2, col = "grey40") +
    geom_vline(xintercept = .start + lubridate::ddays(196), lty = 2, col = "grey40") +
    geom_vline(xintercept = .start + lubridate::ddays(226), lty = 2, col = "grey40") +
    geom_vline(xintercept = .start + lubridate::ddays(251), lty = 2, col = "grey40")
  
  ggsave(paste0("figs-ms/prevalence-", region_names[region], ".pdf"), width = 7, height = 3)
  
}


fcountsp <- reshape2::melt(mprov$post$y_rep) %>%
  #dplyr::rename(day = Var2)
  dplyr::rename(day = Var3) 

ctempp <- fcountsp %>%
  dplyr::filter(day <= ncol(daily_diffs)) %>%
  group_by(day) %>%
  summarise(
    lwr = quantile(value, probs = 0.05),
    lwr2 = quantile(value, probs = 0.25),
    upr = quantile(value, probs = 0.95),
    upr2 = quantile(value, probs = 0.75),
    med = median(value)
  ) %>% mutate(day = .start + lubridate::ddays(day))

flambdasp <- reshape2::melt(mprov$post$lambda_d) %>%
  #dplyr::rename(day = Var2)
  dplyr::rename(day = Var3)

ltempp <- flambdasp %>%
  dplyr::filter(day <= ncol(daily_diffs)) %>%
  group_by(day) %>%
  summarise(
    med = median(value)
  ) %>% mutate(day = .start + lubridate::ddays(day)) %>% mutate(actual = colSums(daily_diffs))

g_countsp <- ggplot(ctempp, aes(x = day, y = med, ymin = lwr, ymax = upr, colour = .hist_blue, fill = .hist_blue)) +
  ylab("Modelled and actual counts") +
  ggtitle("Province-wide") + 
  coord_cartesian(expand = FALSE, xlim = c(.start, .start + lubridate::ddays(mprov$last_day_obs)), ylim = c(0, max(ctempp$upr, colSums(daily_diffs)))) +
  xlab("") + 
  geom_line(data = ltempp, aes(x = day, y = med), alpha = 1, lwd = 1, inherit.aes = FALSE, colour=.hist_blue) +
  geom_ribbon(data = ctempp, alpha = 0.2, fill = .hist_blue, colour = NA, mapping = aes(ymin = lwr, ymax = upr)) +
  geom_ribbon(data = ctempp, alpha = 0.2, mapping = aes(ymin = lwr2, ymax = upr2), fill = .hist_blue, colour = NA) +
  geom_point(data = ltempp, aes(x = day, y = actual), size=1, inherit.aes = FALSE, shape = 1)


library(tidyquant)
## Counts plot
for (region in 1:5) {
  
  fcounts <- reshape2::melt(m7$post$y_rep) %>%
    dplyr::rename(day = Var3) %>%
    dplyr::filter(Var2 %in% region)
  
  ctemp <- fcounts %>%
    dplyr::filter(day <= ncol(daily_diffs)) %>%
    group_by(day) %>%
    summarise(
      lwr = quantile(value, probs = 0.05),
      lwr2 = quantile(value, probs = 0.25),
      upr = quantile(value, probs = 0.95),
      upr2 = quantile(value, probs = 0.75),
      med = median(value)
    ) %>% mutate(day = .start + lubridate::ddays(day))
  
  flambdas <- reshape2::melt(m7$post$lambda_d) %>%
    dplyr::rename(day = Var3)  %>%
    dplyr::filter(Var2 %in% region)
  
  ltemp <- flambdas %>%
    dplyr::filter(day <= ncol(daily_diffs)) %>%
    group_by(day) %>%
    summarise(
      med = median(value)
    ) %>% mutate(day = .start + lubridate::ddays(day)) %>% mutate(actual = daily_diffs[region,])
  
  g_counts <- ggplot(ctemp, aes(x = day, y = med, ymin = lwr, ymax = upr, colour = .hist_blue, fill = .hist_blue)) +
    ylab("Modelled and actual counts") +
    ggtitle(region_names[region]) + 
    coord_cartesian(expand = FALSE, clip = "off", xlim = c(.start, .start + lubridate::ddays(m7$last_day_obs+1)), ylim = c(0, max(ctemp$upr, daily_diffs[region,]))) +
    scale_y_sqrt() + 
    xlab("") + 
    geom_line(data = ltemp, aes(x = day, y = med), alpha = 1, lwd = 1, inherit.aes = FALSE, colour=.hist_blue) +
    geom_ribbon(data = ctemp, alpha = 0.2, fill = .hist_blue, colour = NA, mapping = aes(ymin = lwr, ymax = upr)) +
    geom_ribbon(data = ctemp, alpha = 0.2, mapping = aes(ymin = lwr2, ymax = upr2), fill = .hist_blue, colour = NA) +
    geom_point(data = ltemp, aes(x = day, y = actual), size=1, inherit.aes = FALSE, shape = 1) +
    geom_ma(data = ltemp, aes(x = day, y = actual), inherit.aes = FALSE, ma_fun = SMA, n = 7, lwd = 1, colour='dark red', lty = 1, alpha = 0.5)
    
  
  ggsave(paste0("figs-ms/counts-", region_names[region], ".pdf"), width = 7, height = 3)
  
}


f2p <- mprov$post$f2
f3p <- mprov$post$f3
f4p <- mprov$post$f4
f5p <- mprov$post$f5
f6p <- mprov$post$f6
f7p <- mprov$post$f7
sampFrac_1p <- mprov$post$sampFrac_1
sampFrac_2p <- mprov$post$sampFrac_2
sampFrac_3p <- mprov$post$sampFrac_3
sampFrac_4p <- mprov$post$sampFrac_4
R0p <- tibble(f=mprov$post$R0, para = 1)
R0 <- tibble(f=m7$post$R0, para = 1)

ggplot(R0, aes(x=f)) + geom_density(alpha = 0.7, fill = .hist_blue, colour = NA, adjust = 1.25) +
  geom_density(data = R0p, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
  geom_histogram(aes(y = ..density..),  alpha = 0.3, colour="white", lwd=0.2) +
  coord_cartesian(xlim = c(2.5,3.5), expand = FALSE) +
  xlab("") +   
  scale_x_continuous(breaks = seq(2.6, 3.4, 0.2)) + 
  ylab("Density") + facet_wrap(~para, labeller = label_bquote(R[0][b]))

ggsave(paste0("figs-ms/R0.pdf"), width = 3.2, height = 2.7)

library(latex2exp)
for (region in 1:5) {
  f2 <- m7$post$f2[,region]
  f3 <- m7$post$f3[,region]
  f4 <- m7$post$f4[,region]
  f5 <- m7$post$f5[,region]
  f6 <- m7$post$f6[,region]
  f7 <- m7$post$f7[,region]
  
  com1 <-tibble(f = c(f2,f3,f4,f5,f6,f7), para = rep(c(2,3,4,5,6,7), each=length(f2)))
  com2 <- tibble(f = c(f2p,f3p,f4p,f5p,f6p,f7p), para = rep(c(2,3,4,5,6,7), each=length(f2p)))

  ggplot(com1, aes(x=f)) + geom_density(alpha = 0.7, fill = .hist_blue, colour = NA, adjust = 1.25) +
    geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    geom_histogram(aes(y = ..density..),  alpha = 0.3, colour="white", lwd=0.2) +
    coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab("Density") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(f[ .(para)])) +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/f-", region_names[region], ".pdf"), width = 7, height = 2.7)
  
  sampFrac_1<-m7$post$sampFrac_1[,region]
  sampFrac_2<-m7$post$sampFrac_2[,region]
  sampFrac_3<-m7$post$sampFrac_3[,region]
  sampFrac_4<-m7$post$sampFrac_4[,region]
  com1 <-tibble(f = c(sampFrac_1,sampFrac_2,sampFrac_3,sampFrac_4), R0 = rep(R0$f,4), para = rep(c(1,2,3,4), each=length(sampFrac_1)))
  com2 <- tibble(f = c(sampFrac_1p,sampFrac_2p,sampFrac_3p,sampFrac_4p), para = rep(c(1,2,3,4), each=length(sampFrac_1p)))
  
  ggplot(com1, aes(x=f)) + geom_density(alpha = 1, fill = .hist_blue, colour = NA, adjust = 1.25) + 
    geom_histogram(aes(y = ..density..), fill = .hist_blue, alpha = 0.3, colour="white", lwd=0) + 
    geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab("Density") +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(psi[ .(para)]))  +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/psi_r-", region_names[region], ".pdf"), width = 3.4, height = 2.7)
  
  ggplot(com1, aes(x=f, y=R0)) + geom_point(size = 0.1) + 
    #geom_density(data = com2, colour="NA", fill="dark red", alpha = 0.5, adjust = 1.25) + 
    #coord_cartesian(xlim = c(0,1), expand = FALSE) +
    xlab("") +   
    ylab(TeX("$R_{0b}$")) +
    scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) + facet_wrap(~para, labeller = label_bquote(psi[ .(para)]))  +
    ggtitle(region_names[region]) 
  
  ggsave(paste0("figs-ms/R0-vs-psi_r-", region_names[region], ".pdf"), width = 5.7, height = 2.7)
  
}


print(m7$fit, pars = c("R0", "f2", "f3", "f4", "f5", "f6", "f7", "phi", "lp__", "sampFrac_1", "sampFrac_2", "sampFrac_3", "sampFrac_4"))
#print(mprov$fit, pars = c("R0", "f2", "f3", "f4", "f5", "f6", "f7", "phi", "lp__", "sampFrac_1", "sampFrac_2", "sampFrac_3", "sampFrac_4"))

fit_array<-as.array(m7$fit)

# Get true params if simulation study
source("analysis/simulation-params.R")

my_hist_recover <- function (x, true, facet_args = list(), binwidth = NULL, 
          breaks = NULL) 
{
  x <- bayesplot:::merge_chains(bayesplot:::prepare_mcmc_array(x))
  stopifnot(is.numeric(true), ncol(x) == length(true))
  vline_data <- data.frame(Parameter = colnames(x), True = true)
  lower_data <- data.frame(Parameter = colnames(x), lower = apply(x,2, function(y) quantile(y, 0.05)))
  row.names(lower_data) <- NULL
  upper_data <- data.frame(Parameter = colnames(x), upper = apply(x,2, function(y) quantile(y, 0.95)))
  row.names(upper_data) <- NULL  
  #show(vline_data)
  #show(lower_data)
  hist_data <- bayesplot:::melt_mcmc.matrix(x)[, -1]
  facet_args[["facets"]] <- "Parameter"
  facet_args[["scales"]] <- "free"
  ggplot() + geom_histogram(aes_(x = ~Value, fill = "Estimated"), 
                            data = hist_data, color = bayesplot:::get_color("lh"), size = 0.25, 
                            binwidth = binwidth, breaks = breaks) + 
                            geom_vline(aes_(xintercept = ~True, color = "True"), data = vline_data, size = 1) +
                            geom_vline(aes_(xintercept = ~lower, color = 'red'), data = lower_data, size = 0.5) +
    geom_vline(aes_(xintercept = ~upper, color = 'red'), data = upper_data, size = 0.5) +
    do.call("facet_wrap", facet_args) + scale_fill_manual("", values = c(bayesplot:::get_color("l"), 'blue')) + 
    scale_color_manual("", values = c('red', 'blue')) + bayesplot:::dont_expand_y_axis() + 
    bayesplot::bayesplot_theme_get() + bayesplot:::reduce_legend_spacing(0.25) + 
    bayesplot:::xaxis_title(FALSE) + bayesplot:::yaxis_text(FALSE) + bayesplot:::yaxis_ticks(FALSE) + 
    bayesplot:::yaxis_title(FALSE) + guides(fill="none", color="none") 
}


trace1 <- fit_array[,,c("sampFrac_1[2]", "sampFrac_1[5]", "sampFrac_1[4]", "sampFrac_1[1]", "sampFrac_1[3]",
                        "sampFrac_2[2]", "sampFrac_2[5]", "sampFrac_2[4]", "sampFrac_2[1]", "sampFrac_2[3]",
                        "sampFrac_3[2]", "sampFrac_3[5]", "sampFrac_3[4]", "sampFrac_3[1]", "sampFrac_3[3]",
                        "sampFrac_4[2]", "sampFrac_4[5]", "sampFrac_4[4]", "sampFrac_4[1]", "sampFrac_4[3]")]
dimnames(trace1)[[3]] <- c("ψ1[Coastal]", "ψ1[Fraser]", "ψ1[Interior]", "ψ1[Island]", "ψ1[Northern]",
                           "ψ2[Coastal]", "ψ2[Fraser]", "ψ2[Interior]", "ψ2[Island]", "ψ2[Northern]",
                           "ψ3[Coastal]", "ψ3[Fraser]", "ψ3[Interior]", "ψ3[Island]", "ψ3[Northern]",
                           "ψ4[Coastal]", "ψ4[Fraser]", "ψ4[Interior]", "ψ4[Island]", "ψ4[Northern]")
g <- bayesplot::mcmc_trace(trace1,facet_args = list(dir="v", nrow=5, ncol=4))
ggsave(paste0("figs-ms/traceplot-multiregion-p1.pdf"), width = 8.5, height = 11, device = cairo_pdf)

if (!realdata) {
  g <- my_hist_recover(trace1, c(sampFrac_1[2], sampFrac_1[5], sampFrac_1[4], sampFrac_1[1], sampFrac_1[3],
                                    sampFrac_2[2], sampFrac_2[5], sampFrac_2[4], sampFrac_2[1], sampFrac_2[3],
                                    sampFrac_3[2], sampFrac_3[5], sampFrac_3[4], sampFrac_3[1], sampFrac_3[3],
                                    sampFrac_4[2], sampFrac_4[5], sampFrac_4[4], sampFrac_4[1], sampFrac_4[3]), facet_args = list(dir="v", nrow=5, ncol=4))
  ggsave(paste0("figs-ms/simstudy-multiregion-p1.pdf"), width = 8.5, height = 11, device = cairo_pdf)
}
  
trace1 <- fit_array[,,c("f2[2]", "f2[5]", "f2[4]", "f2[1]", "f2[3]",
                        "f3[2]", "f3[5]", "f3[4]", "f3[1]", "f3[3]",
                        "f4[2]", "f4[5]", "f4[4]", "f4[1]", "f4[3]",
                        "f5[2]", "f5[5]", "f5[4]", "f5[1]", "f5[3]")]
dimnames(trace1)[[3]] <- c("f2[Coastal]", "f2[Fraser]", "f2[Interior]", "f2[Island]", "f2[Northern]",
                           "f3[Coastal]", "f3[Fraser]", "f3[Interior]", "f3[Island]", "f3[Northern]",
                          "f4[Coastal]", "f4[Fraser]", "f4[Interior]", "f4[Island]", "f4[Northern]",
                          "f5[Coastal]", "f5[Fraser]", "f5[Interior]", "f5[Island]", "f5[Northern]")
g <- bayesplot::mcmc_trace(trace1,facet_args = list(dir="v", nrow=5, ncol=4))
ggsave(paste0("figs-ms/traceplot-multiregion-p2.pdf"), width = 8.5, height = 11, device = cairo_pdf)
if (!realdata) {
  g <- my_hist_recover(trace1, c(f2[2], f2[5], f2[4], f2[1], f2[3],
                                              f3[2], f3[5], f3[4], f3[1], f3[3],
                                              f4[2], f4[5], f4[4], f4[1], f4[3],
                                              f5[2], f5[5], f5[4], f5[1], f5[3]), facet_args = list(dir="v", nrow=5, ncol=4))
  ggsave(paste0("figs-ms/simstudy-multiregion-p2.pdf"), width = 8.5, height = 11, device = cairo_pdf)
}


trace1 <- fit_array[,,c("f6[2]", "f6[5]", "f6[4]", "f6[1]", "f6[3]",
                        "f7[2]", "f7[5]", "f7[4]", "f7[1]", "f7[3]",
                        "phi[2]", "phi[5]", "phi[4]", "phi[1]", "phi[3]","R0")]
dimnames(trace1)[[3]] <- c("f6[Coastal]", "f6[Fraser]", "f6[Interior]", "f6[Island]", "f6[Northern]",
                          "f7[Coastal]", "f7[Fraser]", "f7[Interior]", "f7[Island]", "f7[Northern]",
                          "φ[Coastal]", "φ[Fraser]", "φ[Interior]", "φ[Island]", "φ[Northern]","R0b")

g <- bayesplot::mcmc_trace(trace1,facet_args = list(dir="v", nrow=5, ncol=4))
ggsave(paste0("figs-ms/traceplot-multiregion-p3.pdf"), width = 8.5, height = 11, device = cairo_pdf)

if (!realdata) {
  dimnames(trace1)[[3]][16] <- "ω"
  p3labeller <- function(variable,value){
    ifelse(value=='ω',"R0b",value)
  }
  g <- my_hist_recover(trace1, c(f6[2], f6[5], f6[4], f6[1], f6[3],
                                              f7[2], f7[5], f7[4], f7[1], f7[3],
                                              theta[2], theta[5], theta[4], theta[1], theta[3],R0b), facet_args = list(labeller = p3labeller, dir="v", nrow=5, ncol=4))
  ggsave(paste0("figs-ms/simstudy-multiregion-p3.pdf"), width = 8.5, height = 11, device = cairo_pdf)
}


pfit_array<-as.array(mprov$fit)
#tracep <- pfit_array[,,c("sampFrac_1", "sampFrac_2", "sampFrac_3", "sampFrac_4",  "f2", "f3", "f4", "f5", "f6", "f7", "phi[1]", "R0")]
tracep <- pfit_array[,,c("sampFrac_1[1]", "sampFrac_2[1]", "sampFrac_3[1]", "sampFrac_4[1]",  "f2[1]", "f3[1]", "f4[1]", "f5[1]", "f6[1]", "f7[1]", "phi[1]", "R0")]
#tracep <- pfit_array[,,c("R0", "f2", "f3", "f4", "f5", "f6", "phi[1]", "sampFrac_1", "sampFrac_2", "sampFrac_3", "sampFrac_4")]
dimnames(tracep)[[3]] <- c( "ψ1", "ψ2", "ψ3", "ψ4", "f2", "f3", "f4", "f5", "f6", "f7", "φ", "R0b")
g <- bayesplot::mcmc_trace(tracep, facet_args = list(dir="v", nrow=3, ncol=4))
ggsave(paste0("figs-ms/traceplot-provincewide.pdf"), width = 11, height = 8.5, device = cairo_pdf)

