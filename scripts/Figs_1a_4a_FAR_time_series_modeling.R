## Cod FAR
## model the FAR time series to account for model error

library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
obs <- read.csv(".data/ERSST FAR.csv", row.names=1)
mod <- read.csv(".data/CMIP FAR.csv", row.names=1)

obs$year_fac <- as.factor(obs$year)
obs$model_fac <- as.factor(obs$model)

mod$year_fac <- as.factor(mod$year)
mod$model_fac <- as.factor(mod$source)

# change FAR = 1 to Far = 0.9999 to allow beta distribution
change <- obs$FAR == 1 
obs$FAR[change] <- 0.9999

change <- mod$FAR == 1 
mod$FAR[change] <- 0.9999

## Check distribution --------------------------
hist(obs$FAR, breaks = 50)

## brms: setup ---------------------------------------------
## This is the observed FAR time series for Fig. 1a
## Define model formulas
far_formula_fixef <-  bf(FAR ~ year_fac + (1 | model_fac))

## fit: brms --------------------------------------

## observed time series
obs_far_fixef <- brm(far_formula_fixef,
                     data = obs,
                     family = Beta(),
                     cores = 4, chains = 4, iter = 6000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 15))
obs_far_fixef  <- add_criterion(obs_far_fixef, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(obs_far_fixef, file = "output/obs_far_fixef.rds")

obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")
check_hmc_diagnostics(obs_far_fixef$fit)
neff_lowest(obs_far_fixef$fit)
rhat_highest(obs_far_fixef$fit)
summary(obs_far_fixef)
bayes_R2(obs_far_fixef)
y <- obs$FAR
yrep_obs_far_fixef  <- fitted(obs_far_fixef, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_obs_far_fixef[sample(nrow(yrep_obs_far_fixef), 25), ]) +
  xlim(0, 500) +
  ggtitle("obs_far_fixef")

## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.95)[3:4] <- c("ymin.95", "ymax.95")

## 90% CI
ce1s_2 <- conditional_effects(obs_far_fixef, probs = c(0.05, 0.95))
obs.90 <- ce1s_2$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.90)[3:4] <- c("ymin.90", "ymax.90")

## 80% CI
ce1s_3 <- conditional_effects(obs_far_fixef, probs = c(0.1, 0.9))
obs.80 <- ce1s_3$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(obs.80)[3:4] <- c("ymin.80", "ymax.80")


pred.obs <- left_join(obs.95, obs.90)
pred.obs <- left_join(pred.obs, obs.80)
pred.obs$year <- as.numeric(as.character(pred.obs$year_fac))

theme_set(theme_bw())

g1 <- ggplot(pred.obs) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  geom_hline(yintercept = 0.95, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  scale_x_continuous(breaks=seq(1960, 2020, 10)) 

print(g1)

ggsave("./figs/year_predicted_effect_obs_far.png", width = 4.5, height = 2)

## modeled time series of CMIP5 projections for Fig. 4a
## model outputs
mod_far_fixef <- brm(far_formula_fixef,
                     data = mod,
                     family = Beta(),
                     cores = 4, chains = 4, iter = 6000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.99, max_treedepth = 15))
mod_far_fixef  <- add_criterion(mod_far_fixef, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(mod_far_fixef, file = "output/mod_far_fixef.rds")

mod_far_fixef <- readRDS("./output/mod_far_fixef.rds")
check_hmc_diagnostics(mod_far_fixef$fit)
neff_lowest(mod_far_fixef$fit)
rhat_highest(mod_far_fixef$fit)
summary(mod_far_fixef)
bayes_R2(mod_far_fixef)
y <- mod$FAR
yrep_mod_far_fixef  <- fitted(mod_far_fixef, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_mod_far_fixef[sample(nrow(yrep_mod_far_fixef), 25), ]) +
  xlim(0, 500) +
  ggtitle("mod_far_fixef")

## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(mod_far_fixef, probs = c(0.025, 0.975))
mod.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")
mod.95$year <- as.numeric(as.character(mod.95$year_fac))

theme_set(theme_bw())

g1 <- ggplot(mod.95) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90", alpha=0.8) +
  # geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  # geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  geom_hline(yintercept = 0.9, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  scale_x_continuous(breaks=seq(1980, 2040, 10)) 

print(g1)

ggsave("./figs/year_predicted_effect_mod_far.png", width = 4.5, height = 2)

## and a version of the CMIP projections fit w/ 80, 90, 95% CI

## 95% CI
ce1s_1 <- conditional_effects(mod_far_fixef, probs = c(0.025, 0.975))
mod.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")

## 90% CI
ce1s_2 <- conditional_effects(mod_far_fixef, probs = c(0.05, 0.95))
mod.90 <- ce1s_2$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(mod.90)[3:4] <- c("ymin.90", "ymax.90")

## 80% CI
ce1s_3 <- conditional_effects(mod_far_fixef, probs = c(0.1, 0.9))
mod.80 <- ce1s_3$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(mod.80)[3:4] <- c("ymin.80", "ymax.80")


pred.mod <- left_join(mod.95, mod.90)
pred.mod <- left_join(pred.mod, mod.80)
pred.mod$year <- as.numeric(as.character(pred.mod$year_fac))

pred.mod <- pred.mod %>%
  filter(year >= 2006) # limit to RCP8.5 projections

theme_set(theme_bw())

CMIP.FAR <- ggplot(pred.mod) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  # geom_hline(yintercept = 0.90, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("Fraction of attributable risk") +
  scale_x_continuous(breaks=seq(1990, 2040, 10)) 

print(CMIP.FAR)

ggsave("./figs/year_predicted_effect_mod_far_80_90_95_CI.png", width = 4, height = 3)

