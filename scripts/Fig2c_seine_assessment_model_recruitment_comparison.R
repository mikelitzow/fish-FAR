## Model up seine recruitment estimates for each year
## and compare with stock assessment model estimated recruitment
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

theme_set(theme_bw())

## Read in data --------------------------------------------
cod.data <- read.csv("data/cpue.data.csv")
cod.data$cod <- cod.data$cod.age.0
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
                   seed = 123,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(recr_1_zinb, file = "output/recr_1_zinb.rds")
recr_1_zinb  <- add_criterion(recr_1_zinb, c("loo", "bayes_R2"), moment_match = TRUE)

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
                   seed = 12345,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))
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


## Model comparison ---------------------------------------
recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")

loo(recr_1_zinb, recr_2_zinb)

# model 2 best


## Predicted effects ---------------------------------------

## year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

# ce1s_1 is far better!
print(ce1s_1)

# make a plot for the SI
plot.dat <- ce1s_1$year_fac

ggplot(plot.dat, aes(year_fac, estimate__)) +
  geom_point() +
  geom_errorbar(aes(ymax=upper__, ymin=lower__), width = 0.4) +
  coord_trans(y = "pseudo_log") +
  scale_y_continuous(breaks = c(0,1,5,10,50,100,200,500)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Fish / set")

ggsave("./figs/SI_annual_recruitment_estimates_recr_2_zinb.png", width = 6, height = 3.5)

## save the time series to share with collaborators
save.dat <- plot.dat %>%
  select(7:11)

names(save.dat) <- c("year", "cpue", "SE", "LCI", "UCI")
write.csv(save.dat, "./output/cod_seine_annual_cpue_estimates.csv", row.names = F)

# load assessment time series
recr <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv", row.names = 1)

# limiting to 2016 and earlier 
plot <- data.frame(year=2006:2016,
                   ln_assessment_model_R=log(recr$codR0.2020[row.names(recr) %in% 2006:2016]),
                   ln_seine_cpue=log(ce1s_1$year_fac$estimate[1:length(2006:2016)]))

cor(plot) # r = 0.85

ggplot(plot, aes(ln_seine_cpue, ln_assessment_model_R)) +
  geom_text(aes(label=year)) +
  ggtitle("2006-2016") 

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
                 seed = 1234,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.999, max_treedepth = 10))
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
                 seed = 1234,
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
