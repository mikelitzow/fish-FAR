## Model up seine recruitment estimates for each year
## and compare with stock assessment model esimated recruitment

library(ggplot2)
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
recr <- read.csv("data/stock_assessment_recruitment.csv", row.names = 1)

plot <- data.frame(year=2006:2017,
                   ln_assessment_model_R=log(recr$recruits[row.names(recr) %in% 2006:2017]),
                   ln_seine_cpue=log(ce1s_1$year_fac$estimate[1:length(2006:2017)]))

cor(plot) # r = 0.678

ggplot(plot, aes(ln_seine_cpue, ln_assessment_model_R)) +
  geom_text(label=year) +
  ggtitle("2006-2017") +
  theme_bw()

ggsave("./figs/seine_assessment_model_recruitment_comparison.png", width = 5, height = 4)
