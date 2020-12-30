## correlation of egg/larval temperature anomalies with
## ERSST anomalies for area covered by Walsh et al. 2018 and
## Walsh et al. Fraction of Attributable Risk estimates

library(tidyverse)
library(brms)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# compare ERSST anomalies with GODAS cod/larvae anomalies
ersst <- read.csv("climate/GOA ERSST anomalies matching Walsh et al area.csv")
godas <- read.csv("climate/godas anomalies.csv")

godas <- left_join(godas, ersst)
ggplot(godas, aes(anomaly, mean.anom)) +
  geom_point()

cor(godas$anomaly, godas$mean.anom) # 0.737

# this is the estimate of FAR across all five models for each year
obs_far <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

far <- left_join(godas, obs)

ggplot(far, aes(estimate__, mean.anom)) +
  geom_point()

cor(far$estimate__, far$mean.anom) # r = 0.761