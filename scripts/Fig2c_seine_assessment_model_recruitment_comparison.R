## Model up seine recruitment estimates for each year
## and compare with stock assessment model esimated recruitment
## 
## Used to estimate stock assessment model estimates for years not supported by data
## in the model (2017-2020) for Fig. 2c

## Also includes extensions for recruitment prediction offshoot idea
library(tidyverse)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
cod.data <- read.csv("data/cpue.data.csv", row.names = 1)
cod.data$bay_fac <- as.factor(cod.data$bay)
cod.data$year_fac <- as.factor(cod.data$year)
cod.data$site_fac <- as.factor(cod.data$site)
cod.data$bay_site_fac <- as.factor(paste0(cod.data$bay, "_", cod.data$site))
cod.data$present <- ifelse(cod.data$cod > 0, 1, 0)
cod.data$date <- as.Date(cod.data$julian,
                         origin = paste0(cod.data$year, "-01-01"))

## brms: setup ---------------------------------------------

## Define model formulas
recr_1_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac) + (year_fac),
                    zi ~ s(julian, k = 3) + (1 | bay_fac) + (year_fac))

recr_2_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac))

## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Show default priors
get_prior(recr_1_formula, cod.data, family = zinb)


## Set priors
priors_zinb <- c(set_prior("normal(0, 3)", class = "b"),
                    set_prior("normal(0, 3)", class = "Intercept"),
                    set_prior("student_t(3, 0, 3)", class = "sd"),
                    set_prior("student_t(3, 0, 3)", class = "sds"),
                    set_prior("gamma(0.01, 0.01)", class = "shape"),
                    set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                    set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                    set_prior("student_t(3, 0, 3)", class = "sd", dpar = "zi"),
                    set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))


## fit: zero-inflated --------------------------------------
recr_1_zinb <- brm(recr_1_formula,
                    data = cod.data,
                    prior = priors_zinb,
                    family = zinb,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
recr_1_zinb  <- add_criterion(recr_1_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(recr_1_zinb, file = "output/recr_1_zinb.rds")

recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")
check_hmc_diagnostics(recr_1_zinb$fit)
neff_lowest(recr_1_zinb$fit)
rhat_highest(recr_1_zinb$fit)
summary(recr_1_zinb)
bayes_R2(recr_1_zinb)
plot(recr_1_zinb$criteria$loo, "k")
plot(conditional_smooths(recr_1_zinb), ask = FALSE)
plot(recr_1_zinb, ask = FALSE)
y <- cod.data$cod
yrep_recr_1_zinb  <- fitted(recr_1_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_recr_1_zinb[sample(nrow(yrep_recr_1_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("recr_1_zinb")
pdf("./figs/trace_recr_1_zinb.pdf", width = 6, height = 4)
trace_plot(recr_1_zinb$fit)
dev.off()

recr_2_zinb <- brm(recr_2_formula,
                   data = cod.data,
                   prior = priors_zinb,
                   family = zinb,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
recr_2_zinb  <- add_criterion(recr_2_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(recr_2_zinb, file = "output/recr_2_zinb.rds")

recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")
check_hmc_diagnostics(recr_2_zinb$fit)
neff_lowest(recr_2_zinb$fit)
rhat_highest(recr_2_zinb$fit)
summary(recr_2_zinb)
bayes_R2(recr_2_zinb)
plot(recr_2_zinb$criteria$loo, "k")
plot(conditional_smooths(recr_2_zinb), ask = FALSE)
y <- cod.data$cod
yrep_recr_2_zinb  <- fitted(recr_2_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_recr_2_zinb[sample(nrow(yrep_recr_2_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("recr_2_zinb")
pdf("./figs/trace_recr_2_zinb.pdf", width = 6, height = 4)
trace_plot(recr_2_zinb$fit)
dev.off()

## fit: only 2006-2017-----------------------------------
recr_1_zinb_06_17 <- brm(recr_1_formula,
                   data = cod.data[cod.data$year %in% 2006:2017,],
                   prior = priors_zinb,
                   family = zinb,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
recr_1_zinb_06_17  <- add_criterion(recr_1_zinb_06_17, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(recr_1_zinb_06_17, file = "output/recr_1_zinb_06_17.rds")

recr_1_zinb_06_17 <- readRDS("./output/recr_1_zinb_06_17.rds")
check_hmc_diagnostics(recr_1_zinb_06_17$fit)
neff_lowest(recr_1_zinb_06_17$fit)
rhat_highest(recr_1_zinb_06_17$fit)
summary(recr_1_zinb_06_17)
bayes_R2(recr_1_zinb_06_17)
plot(recr_1_zinb_06_17$criteria$loo, "k")
plot(conditional_smooths(recr_1_zinb_06_17), ask = FALSE)
plot(recr_1_zinb_06_17, ask = FALSE)
y <- cod.data$cod
yrep_recr_1_zinb_06_17  <- fitted(recr_1_zinb_06_17, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_recr_1_zinb_06_17[sample(nrow(yrep_recr_1_zinb_06_17), 25), ]) +
  xlim(0, 500) +
  ggtitle("recr_1_zinb_06_17")
pdf("./figs/trace_recr_1_zinb_06_17.pdf", width = 6, height = 4)
trace_plot(recr_1_zinb_06_17$fit)
dev.off()

recr_2_zinb_06_17 <- brm(recr_2_formula,
                   data = cod.data[cod.data$year %in% 2006:2017,],
                   prior = priors_zinb,
                   family = zinb,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
recr_2_zinb_06_17  <- add_criterion(recr_2_zinb_06_17, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(recr_2_zinb_06_17, file = "output/recr_2_zinb_06_17.rds")

recr_2_zinb_06_17 <- readRDS("./output/recr_2_zinb_06_17.rds")
check_hmc_diagnostics(recr_2_zinb_06_17$fit)
neff_lowest(recr_2_zinb_06_17$fit)
rhat_highest(recr_2_zinb_06_17$fit)
summary(recr_2_zinb_06_17)
bayes_R2(recr_2_zinb_06_17)
plot(recr_2_zinb_06_17$criteria$loo, "k")
plot(conditional_smooths(recr_2_zinb_06_17), ask = FALSE)
y <- cod.data$cod
yrep_recr_2_zinb_06_17  <- fitted(recr_2_zinb_06_17, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_recr_2_zinb_06_17[sample(nrow(yrep_recr_2_zinb_06_17), 25), ]) +
  xlim(0, 500) +
  ggtitle("recr_2_zinb_06_17")
pdf("./figs/trace_recr_2_zinb_06_17.pdf", width = 6, height = 4)
trace_plot(recr_2_zinb_06_17$fit)
dev.off()


## Model comparison ---------------------------------------
recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")
recr_1_zinb_06_17 <- readRDS("./output/recr_1_zinb_06_17.rds")
recr_2_zinb_06_17 <- readRDS("./output/recr_2_zinb_06_17.rds")

loo(recr_1_zinb, recr_2_zinb)
loo(recr_1_zinb_06_17, recr_2_zinb_06_17)

# model 2 best in both cases


## Predicted effects ---------------------------------------

## year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

# ce1s_1 is far better!
print(ce1s_1)
ggsave("./figs/annual_recruitment_estimates_recr_2_zinb.png", width = 7, height = 4)

# load assessment time series
recr <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv", row.names = 1)

# limiting to 2016 and earlier 
plot <- data.frame(year=2006:2016,
                   ln_assessment_model_R=log(recr$codR0.2020[row.names(recr) %in% 2006:2016]),
                   ln_seine_cpue=log(ce1s_1$year_fac$estimate[1:length(2006:2016)]))

cor(plot) # r = 0.86

ggplot(plot, aes(ln_seine_cpue, ln_assessment_model_R)) +
  geom_text(aes(label=year)) +
  ggtitle("2006-2016") +
  theme_bw()

ggsave("./figs/seine_assessment_model_recruitment_comparison.png", width = 5, height = 4)

## fit a brms model ------------------------------------------

# first, clean up data
dat <- plot

names(dat)[2:3] <- c("model", "seine")

# and add ssb
ssb <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

ssb <- ssb %>%
  select(year, codSSB.2020)

names(ssb)[2] <- "ssb"

dat <- left_join(dat, ssb)

## brms: setup ---------------------------------------------

## Define model formulas
## Limiting knots to 3 to prevent overfitting

codR1_formula <-  bf(model ~ s(ssb, k = 3) + seine)

codR2_formula <-  bf(model ~ seine)


## fit --------------------------------------
codR1_brm <- brm(codR1_formula,
                    data = dat,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
codR1_brm  <- add_criterion(codR1_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR1_brm, file = "output/codR1_brm.rds")

codR1_brm <- readRDS("./output/codR1_brm.rds")
check_hmc_diagnostics(codR1_brm$fit)
neff_lowest(codR1_brm$fit)
rhat_highest(codR1_brm$fit)
summary(codR1_brm)
bayes_R2(codR1_brm)
plot(codR1_brm$criteria$loo, "k")
plot(conditional_smooths(codR1_brm), ask = FALSE)
y <- trend$trend
yrep_codR1_brm  <- fitted(codR1_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR1_brm[sample(nrow(yrep_codR1_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR1_brm")
pdf("./figs/trace_codR1_brm.pdf", width = 6, height = 4)
trace_plot(codR1_brm$fit)
dev.off()


codR2_brm <- brm(codR2_formula,
                     data = dat,
                     cores = 4, chains = 4, iter = 3000,
                     save_pars = save_pars(all = TRUE),
                     control = list(adapt_delta = 0.999, max_treedepth = 10))
codR2_brm  <- add_criterion(codR2_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR2_brm, file = "output/codR2_brm.rds")

codR2_brm <- readRDS("./output/codR2_brm.rds")
check_hmc_diagnostics(codR2_brm$fit)
neff_lowest(codR2_brm$fit)
rhat_highest(codR2_brm$fit)
summary(codR2_brm)
bayes_R2(codR2_brm)
plot(codR2_brm$criteria$loo, "k")
plot(conditional_smooths(codR2_brm), ask = FALSE)
y <- trend$trend
yrep_codR2_brm  <- fitted(codR2_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR2_brm[sample(nrow(yrep_codR2_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR2_brm")
pdf("./figs/trace_codR2_brm.pdf", width = 6, height = 4)
trace_plot(codR2_brm$fit)
dev.off()


## Model selection -----------------------------------------
codR1_brm  <- readRDS("./output/codR1_brm.rds")
codR2_brm  <- readRDS("./output/codR2_brm.rds")

loo(codR1_brm, codR2_brm)

## second round of model-fitting adding seine:FAR interaction

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

dat <- left_join(dat, obs)
names(dat)[6] <- "far"

## Define model formulas

codR3_formula <-  bf(model ~ seine + seine:far)

codR4_formula <-  bf(model ~  seine:far)
## fit --------------------------------------
codR3_brm <- brm(codR3_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
codR3_brm  <- add_criterion(codR3_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR3_brm, file = "output/codR3_brm.rds")

codR3_brm <- readRDS("./output/codR3_brm.rds")
check_hmc_diagnostics(codR3_brm$fit)
neff_lowest(codR3_brm$fit)
rhat_highest(codR3_brm$fit)
summary(codR3_brm)
bayes_R2(codR3_brm)
plot(codR3_brm$criteria$loo, "k")
plot(conditional_smooths(codR3_brm), ask = FALSE)
y <- trend$trend
yrep_codR3_brm  <- fitted(codR3_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR3_brm[sample(nrow(yrep_codR3_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR3_brm")
pdf("./figs/trace_codR3_brm.pdf", width = 6, height = 4)
trace_plot(codR3_brm$fit)
dev.off()


codR4_brm <- brm(codR4_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
codR4_brm  <- add_criterion(codR4_brm, c("loo", "bayes_R4"), moment_match = TRUE)
saveRDS(codR4_brm, file = "output/codR4_brm.rds")

codR4_brm <- readRDS("./output/codR4_brm.rds")
check_hmc_diagnostics(codR4_brm$fit)
neff_lowest(codR4_brm$fit)
rhat_highest(codR4_brm$fit)
summary(codR4_brm)
bayes_R4(codR4_brm)
plot(codR4_brm$criteria$loo, "k")
plot(conditional_smooths(codR4_brm), ask = FALSE)
y <- trend$trend
yrep_codR4_brm  <- fitted(codR4_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR4_brm[sample(nrow(yrep_codR4_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR4_brm")
pdf("./figs/trace_codR4_brm.pdf", width = 6, height = 4)
trace_plot(codR4_brm$fit)
dev.off()

## model selection --------------------------------------
codR2_brm  <- readRDS("./output/codR2_brm.rds")
codR3_brm  <- readRDS("./output/codR3_brm.rds")
codR4_brm  <- readRDS("./output/codR4_brm.rds")

loo(codR2_brm, codR3_brm, codR4_brm)

## plot predicted values codR2_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(codR2_brm, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR2_brm, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR2_brm, effect = "seine", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "ln(seine)", y = "ln(model recruitment)") +
  geom_text(data=dat, aes(seine, model, label = year), size=1.5) +
  theme_bw()

print(g)

ggsave("./figs/codR2_brm_seine.png", width=3, height=2.5, units = 'in')


## plot predicted values codR3_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(codR3_brm, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR3_brm, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR3_brm, effect = "seine", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$seine), amount = 0.01),
                          rep(NA, 100-length(unique(dat$seine))))

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "ln(seine)", y = "ln(model recruitment)") +
  theme_bw()+
  geom_rug(aes(x=rug.anom, y=NULL))
print(g)


ggsave("./figs/codR3_brm_seine.png", width=3, height=2, units = 'in')


## and far
## plot predicted values ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(codR3_brm, effect = "far", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(codR3_brm, effect = "far", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(codR3_brm, effect = "far", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$far
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$far[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$far[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$far[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$far[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$far), amount = 0.01),
                          rep(NA, 100-length(unique(dat$far))))

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "FAR", y = "ln(model recruitment)") +
  theme_bw()+
  geom_rug(aes(x=rug.anom, y=NULL))
print(g)

ggsave("./figs/codR3_brm_far.png", width=3, height=2, units = 'in')

## try plotting the seine-far interaction ------------------------------------------
# plot as Z-scores
dat$z.score <- as.vector(scale(dat$far))

hist(dat$z.score)

dat$FAR <- if_else(dat$z.score < -0.5, "< -0.5 SD", 
                   if_else(dat$z.score > 0.5, "> 0.5 SD", "Mean"))

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())

ggplot(dat, aes(seine, model, color=FAR)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  scale_color_manual(values=cb[c(2,4,6)])
## not enough data!
ggsave("./figs/cod_seine_FAR_R.png", width=4.5, height=2.5, units = 'in')


## brm: seine and far - full interactive and main term ------------

codR5_formula <-  bf(model ~  seine*far)

codR5_brm <- brm(codR5_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
codR5_brm  <- add_criterion(codR5_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR5_brm, file = "output/codR5_brm.rds")

codR5_brm <- readRDS("./output/codR5_brm.rds")
check_hmc_diagnostics(codR5_brm$fit)
neff_lowest(codR5_brm$fit)
rhat_highest(codR5_brm$fit)
summary(codR5_brm)
bayes_R2(codR5_brm)
plot(codR5_brm$criteria$loo, "k")
plot(conditional_smooths(codR5_brm), ask = FALSE)
y <- trend$trend
yrep_codR5_brm  <- fitted(codR5_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_codR5_brm[sample(nrow(yrep_codR5_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("codR5_brm")
pdf("./figs/trace_codR5_brm.pdf", width = 6, height = 4)
trace_plot(codR5_brm$fit)
dev.off()

## model selection --------------------------------------
codR2_brm  <- readRDS("./output/codR2_brm.rds")
codR3_brm  <- readRDS("./output/codR3_brm.rds")
codR4_brm  <- readRDS("./output/codR4_brm.rds")
codR5_brm  <- readRDS("./output/codR5_brm.rds")

loo(codR2_brm, codR3_brm, codR4_brm, codR5_brm)