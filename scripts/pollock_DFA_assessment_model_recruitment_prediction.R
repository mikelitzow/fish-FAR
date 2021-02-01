
## Compare DFA trend on three field time series with pollock stock assessment model esimated recruitment

library(tidyverse)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

# load DFA trend
dat <- read.csv("./output/poll_dfa_trend.csv", row.names = 1)

dat <- dat %>%
  select(year, trend)

names(dat)[2] <- "dfa"

# and add ssb
ssb <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

ssb <- ssb %>%
  select(year, poll.SSB.2020, pollR0.2020)

names(ssb)[2:3] <- c("ssb", "model")

dat <- left_join(dat, ssb)

dat$model <- log(dat$model)

## no model recruitment estimates available for 2020 - drop
dat <- dat %>%
  filter(year < 2020)

## Define model formulas
## Limiting knots to 3 to prevent overfitting

pollR1_formula <-  bf(model ~ s(ssb, k = 3) + s(dfa, k = 3))

pollR2_formula <-  bf(model ~ s(dfa, k = 3))


## fit --------------------------------------
pollR1_brm <- brm(pollR1_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
pollR1_brm  <- add_criterion(pollR1_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR1_brm, file = "output/pollR1_brm.rds")

pollR1_brm <- readRDS("./output/pollR1_brm.rds")
check_hmc_diagnostics(pollR1_brm$fit)
neff_lowest(pollR1_brm$fit)
rhat_highest(pollR1_brm$fit)
summary(pollR1_brm)
bayes_R2(pollR1_brm)
plot(pollR1_brm$criteria$loo, "k")
plot(conditional_smooths(pollR1_brm), ask = FALSE)
y <- trend$trend
yrep_pollR1_brm  <- fitted(pollR1_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR1_brm[sample(nrow(yrep_pollR1_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollR1_brm")
pdf("./figs/trace_pollR1_brm.pdf", width = 6, height = 4)
trace_plot(pollR1_brm$fit)
dev.off()


pollR2_brm <- brm(pollR2_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR2_brm  <- add_criterion(pollR2_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR2_brm, file = "output/pollR2_brm.rds")

pollR2_brm <- readRDS("./output/pollR2_brm.rds")
check_hmc_diagnostics(pollR2_brm$fit)
neff_lowest(pollR2_brm$fit)
rhat_highest(pollR2_brm$fit)
summary(pollR2_brm)
bayes_R2(pollR2_brm)
plot(pollR2_brm$criteria$loo, "k")
plot(conditional_smooths(pollR2_brm), ask = FALSE)
y <- trend$trend
yrep_pollR2_brm  <- fitted(pollR2_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR2_brm[sample(nrow(yrep_pollR2_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollR2_brm")
pdf("./figs/trace_pollR2_brm.pdf", width = 6, height = 4)
trace_plot(pollR2_brm$fit)
dev.off()


## Model selection -----------------------------------------
pollR1_brm  <- readRDS("./output/pollR1_brm.rds")
pollR2_brm  <- readRDS("./output/pollR2_brm.rds")

loo(pollR1_brm, pollR2_brm)

## second round of model-fitting adding dfa:FAR interaction

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

dat <- left_join(dat, obs)
names(dat)[6] <- "far"

## Define model formulas

pollR3_formula <-  bf(model ~ s(dfa, k = 3) + dfa:far)

pollR4_formula <-  bf(model ~ dfa:far)


## fit --------------------------------------
pollR3_brm <- brm(pollR3_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
pollR3_brm  <- add_criterion(pollR3_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR3_brm, file = "output/pollR3_brm.rds")

pollR3_brm <- readRDS("./output/pollR3_brm.rds")
check_hmc_diagnostics(pollR3_brm$fit)
neff_lowest(pollR3_brm$fit)
rhat_highest(pollR3_brm$fit)
summary(pollR3_brm)
bayes_R2(pollR3_brm)
plot(pollR3_brm$criteria$loo, "k")
plot(conditional_smooths(pollR3_brm), ask = FALSE)
y <- trend$trend
yrep_pollR3_brm  <- fitted(pollR3_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR3_brm[sample(nrow(yrep_pollR3_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollR3_brm")
pdf("./figs/trace_pollR3_brm.pdf", width = 6, height = 4)
trace_plot(pollR3_brm$fit)
dev.off()


pollR4_brm <- brm(pollR4_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR4_brm  <- add_criterion(pollR4_brm, c("loo", "bayes_R4"), moment_match = TRUE)
saveRDS(pollR4_brm, file = "output/pollR4_brm.rds")

pollR4_brm <- readRDS("./output/pollR4_brm.rds")
check_hmc_diagnostics(pollR4_brm$fit)
neff_lowest(pollR4_brm$fit)
rhat_highest(pollR4_brm$fit)
summary(pollR4_brm)
bayes_R2(pollR4_brm)
plot(pollR4_brm$criteria$loo, "k")
plot(conditional_smooths(pollR4_brm), ask = FALSE)
y <- trend$trend
yrep_pollR4_brm  <- fitted(pollR4_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR4_brm[sample(nrow(yrep_pollR4_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollR4_brm")
pdf("./figs/trace_pollR4_brm.pdf", width = 6, height = 4)
trace_plot(pollR4_brm$fit)
dev.off()

## model selection --------------------------------------
pollR3_brm  <- readRDS("./output/pollR3_brm.rds")
pollR4_brm  <- readRDS("./output/pollR4_brm.rds")

loo(pollR1_brm, pollR2_brm, pollR3_brm, pollR4_brm)

## plot predicted values pollR2_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(pollR2_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR2_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR2_brm, effect = "dfa", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "ln(model recruitment)") +
  geom_text(data=dat, aes(dfa, model, label = year), size=1.5) +
  theme_bw()

print(g)

ggsave("./figs/pollR2_brm_dfa.png", width=3.5, height=2.5, units = 'in')


## try plotting the interaction
library(interactions)

mod <- lm(model ~ dfa:far, data=dat)
summary(mod)

interact_plot(mod, pred = dfa, modx = far, plot.points = T)

## might be better to plot as Z-scores
dat$z.score <- as.vector(scale(dat$far))

hist(dat$z.score)

dat$FAR <- if_else(dat$z.score < -1, "< -1SD", 
                        if_else(dat$z.score > 1, "> 1SD", "Mean"))

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())

ggplot(dat, aes(dfa, model, color=FAR)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values=cb[c(2,4,6)])

ggsave("./figs/poll_DFA_FAR_R.png", width=4.5, height=2.5, units = 'in')

## brm: seine and far - full interactive and main term ------------

pollR5_formula <-  bf(model ~  dfa*far)

pollR5_brm <- brm(pollR5_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
pollR5_brm  <- add_criterion(pollR5_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pollR5_brm, file = "output/pollR5_brm.rds")

pollR5_brm <- readRDS("./output/pollR5_brm.rds")
check_hmc_diagnostics(pollR5_brm$fit)
neff_lowest(pollR5_brm$fit)
rhat_highest(pollR5_brm$fit)
summary(pollR5_brm)
bayes_R2(pollR5_brm)
plot(pollR5_brm$criteria$loo, "k")
plot(conditional_smooths(pollR5_brm), ask = FALSE)
y <- trend$trend
yrep_pollR5_brm  <- fitted(pollR5_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollR5_brm[sample(nrow(yrep_pollR5_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollR5_brm")
pdf("./figs/trace_pollR5_brm.pdf", width = 6, height = 4)
trace_plot(pollR5_brm$fit)
dev.off()

## model selection --------------------------------------
pollR2_brm  <- readRDS("./output/pollR2_brm.rds")
pollR3_brm  <- readRDS("./output/pollR3_brm.rds")
pollR4_brm  <- readRDS("./output/pollR4_brm.rds")
pollR5_brm  <- readRDS("./output/pollR5_brm.rds")

loo(pollR2_brm, pollR3_brm, pollR4_brm, pollR5_brm)
