# Make plots based on multiregion and provincial-wide fits for November prediction
# Needs results saved from "run-multiregion-predNov.R" and "run-provincewide-predNov.R"
# Some plotting code adapted by Samuel Wong from original code by authors of https://github.com/carolinecolijn/distancing-impact-covid19

library(rstan)

Island <- read.csv("data-generated/daily-cases-Island.csv")
Interior <- read.csv("data-generated/daily-cases-Interior.csv")
Coastal <- read.csv("data-generated/daily-cases-Coastal.csv")
Northern <- read.csv("data-generated/daily-cases-Northern.csv")
Fraser <- read.csv("data-generated/daily-cases-Fraser.csv")

daily_diffs <- rbind(Island$cases,Coastal$cases,Northern$cases,Interior$cases,Fraser$cases)
pop_size <- c(843375,1225195,297570,795116,1889225)

load("data-generated/fit-predNov-multiregion-loosepriors-mode.rda")
mprov <- readRDS("data-generated/fit-predNov-provincewide-loosepriors-mode.rds")

region_names <- c("Island", "Coastal", "Northern", "Interior", "Fraser")
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
  dplyr::filter(time < ncol(daily_diffs)+24) %>%
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
    dplyr::filter(time < ncol(daily_diffs)+24) %>%
    group_by(day) %>%
    summarise(
      lwr = quantile(prevalence, probs = 0.05),
      lwr2 = quantile(prevalence, probs = 0.25),
      upr = quantile(prevalence, probs = 0.95),
      upr2 = quantile(prevalence, probs = 0.75),
      med = median(prevalence)
    )

  y_limits <- c(0, max(c(temp$upr, tempp$upr * pop_size[region] / sum(pop_size))) * 1.01)
  g_prev <- ggplot(temp, aes(x = day, y = med, ymin = lwr, ymax = upr, colour = .hist_blue, fill = .hist_blue)) +
    ylab("Modelled prevalence") +
    ggtitle(region_names[region]) + 
    coord_cartesian(expand = FALSE, xlim = c(.start + lubridate::ddays(m7$last_day_obs - 37), .start + lubridate::ddays(m7$last_day_obs + 23)), ylim = y_limits) +
    annotate("rect", xmin = .start + lubridate::ddays(m7$last_day_obs), xmax = .start + lubridate::ddays(m7$last_day_obs + 23), ymin = 0, ymax = y_limits[2], fill = "grey95") +
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
    geom_vline(xintercept = .start + lubridate::ddays(226), lty = 2, col = "grey40") #+
    #geom_vline(xintercept = .start + lubridate::ddays(251), lty = 2, col = "grey40")
  
  ggsave(paste0("figs-ms/prevalence-", region_names[region], "-prediction.pdf"), width = 3, height = 3)
  
}





