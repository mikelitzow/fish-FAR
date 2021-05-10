## Cod analysis - effect of temperature on seine recruitment estimates
## Negative binomial
## Smooths limited to k=3 to avoid over-fitting
## For Fig. 2a in the draft
## note that cod2sg_zinb_k3 is the selected model reported in the ms.

library(ggplot2)
library(dplyr)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


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

# add ssb data
ssb.data <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")
ssb.data <- ssb.data %>%
    select(year, codSSB.2020)
names(ssb.data)[2] <- "ssb"
cod.data <- left_join(cod.data, ssb.data)

# add temperature data
temp.data <- read.csv("data/godas.anomalies.csv")
cod.data$temp.anom <- temp.data$mean.anom[match(as.numeric(as.character(cod.data$year)), temp.data$year)]


## brms: setup ---------------------------------------------

## Define model formulas
cod0_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac),
                     zi ~ s(julian, k = 3) + (1 | bay_fac))

cod0s_formula <-  bf(cod ~ s(julian, k = 3) + (1 | bay_fac/site_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac))

cod1sg_formula <-  bf(cod ~ s(julian, k = 3) + s(temp.anom, k = 3) + (1 | bay_fac/site_fac),
                     zi ~ s(julian, k = 3) + s(temp.anom, k = 3) + (1 | bay_fac/site_fac))

cod2sg_formula <-  bf(cod ~ s(julian, k = 3) + s(temp.anom, k = 3) + s(ssb, k = 3) +
                         (1 | bay_fac/site_fac),
                     zi ~ s(julian, k = 3) + s(temp.anom, k = 3) + s(ssb, k = 3) +
                         (1 | bay_fac/site_fac))


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")


## Show default priors
# get_prior(cod0_formula, cod.data, family = zinb)
# get_prior(cod0_formula_nb, cod.data, family = nb)

## Set priors
priors_zinb_k3 <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 3)", class = "sd"),
                 set_prior("student_t(3, 0, 3)", class = "sds"),
                 set_prior("gamma(0.01, 0.01)", class = "shape"),
                 set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                 set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sd", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))


## fit: zero-inflated --------------------------------------
cod0_zinb_k3 <- brm(cod0_formula,
                 data = cod.data,
                 prior = priors_zinb_k3,
                 family = zinb,
                 seed = 1234,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
cod0_zinb_k3  <- add_criterion(cod0_zinb_k3, c("loo", "bayes_R2"),
                               moment_match = TRUE)
saveRDS(cod0_zinb_k3, file = "output/cod0_zinb_k3.rds")

cod0_zinb_k3 <- readRDS("./output/cod0_zinb_k3.rds")
check_hmc_diagnostics(cod0_zinb_k3$fit)
neff_lowest(cod0_zinb_k3$fit)
rhat_highest(cod0_zinb_k3$fit)
summary(cod0_zinb_k3)
bayes_R2(cod0_zinb_k3)
plot(cod0_zinb_k3$criteria$loo, "k")
plot(conditional_smooths(cod0_zinb_k3), ask = FALSE)
y <- cod.data$cod
yrep_cod0_zinb_k3  <- fitted(cod0_zinb_k3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod0_zinb_k3[sample(nrow(yrep_cod0_zinb_k3), 25), ]) +
    xlim(0, 500) +
    ggtitle("cod0_zinb_k3")
pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
    trace_plot(cod0_zinb_k3$fit)
dev.off()


cod0s_zinb_k3 <- brm(cod0s_formula,
                  data = cod.data,
                  prior = priors_zinb_k3,
                  family = zinb,
                  seed = 1234,
                  cores = 4, chains = 4, iter = 4000,
                  save_pars = save_pars(all = TRUE),
                  control = list(adapt_delta = 0.99, max_treedepth = 10))
cod0s_zinb_k3  <- add_criterion(cod0s_zinb_k3, c("loo", "bayes_R2"),
                                moment_match = TRUE, reloo = FALSE,
                                cores = 4, k_threshold = 0.7)
saveRDS(cod0s_zinb_k3, file = "output/cod0s_zinb_k3.rds")

cod0s_zinb_k3 <- readRDS("./output/cod0s_zinb_k3.rds")
check_hmc_diagnostics(cod0s_zinb_k3$fit)
neff_lowest(cod0s_zinb_k3$fit)
rhat_highest(cod0s_zinb_k3$fit)
summary(cod0s_zinb_k3)
bayes_R2(cod0s_zinb_k3)
plot(cod0s_zinb_k3$criteria$loo, "k")
plot(conditional_smooths(cod0s_zinb_k3), ask = FALSE)
y <- cod.data$cod
yrep_cod0s_zinb_k3  <- fitted(cod0s_zinb_k3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod0s_zinb_k3[sample(nrow(yrep_cod0s_zinb_k3), 25), ]) +
    xlim(0, 500) +
    ggtitle("cod0s_zinb_k3")
pdf("./figs/trace_cod0s_zinb_k3.pdf", width = 6, height = 4)
    trace_plot(cod0s_zinb_k3$fit)
dev.off()


## fit: GODAS anomaly models -------------------------------
## reloo**
cod1sg_zinb_k3 <- brm(cod1sg_formula,
                   data = cod.data,
                   prior = priors_zinb_k3,
                   family = zinb,
                   seed = 1234,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))
cod1sg_zinb_k3  <- add_criterion(cod1sg_zinb_k3, c("loo", "bayes_R2"),
                                 moment_match = TRUE, reloo = TRUE,
                                 cores = 4, k_threshold = 0.7)
saveRDS(cod1sg_zinb_k3, file = "output/cod1sg_zinb_k3.rds")

cod1sg_zinb_k3 <- readRDS("./output/cod1sg_zinb_k3.rds")
check_hmc_diagnostics(cod1sg_zinb_k3$fit)
neff_lowest(cod1sg_zinb_k3$fit)
rhat_highest(cod1sg_zinb_k3$fit)
summary(cod1sg_zinb_k3)
bayes_R2(cod1sg_zinb_k3)
plot(cod1sg_zinb_k3$criteria$loo, "k")
plot(conditional_smooths(cod1sg_zinb_k3), ask = FALSE)
pdf("./figs/trace_cod1sg_zinb_k3.pdf", width = 6, height = 4)
    trace_plot(cod1sg_zinb_k3$fit)
dev.off()


cod2sg_zinb_k3 <- brm(cod2sg_formula,
                   data = cod.data,
                   prior = priors_zinb_k3,
                   family = zinb,
                   seed = 1234,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))
cod2sg_zinb_k3  <- add_criterion(cod2sg_zinb_k3, c("loo", "bayes_R2"),
                                 moment_match = TRUE, reloo = TRUE,
                                 cores = 4, k_threshold = 0.7)
saveRDS(cod2sg_zinb_k3, file = "output/cod2sg_zinb_k3.rds")

cod2sg_zinb_k3 <- readRDS("./output/cod2sg_zinb_k3.rds")
check_hmc_diagnostics(cod2sg_zinb_k3$fit)
neff_lowest(cod2sg_zinb_k3$fit)
rhat_highest(cod2sg_zinb_k3$fit)
summary(cod2sg_zinb_k3)
bayes_R2(cod2sg_zinb_k3)
plot(cod2sg_zinb_k3$criteria$loo, "k")
plot(conditional_smooths(cod2sg_zinb_k3), ask = FALSE)
pdf("./figs/trace_cod2sg_zinb_k3.pdf", width = 6, height = 4)
    trace_plot(cod2sg_zinb_k3$fit)
dev.off()


## Model selection -----------------------------------------
cod0_zinb_k3   <- readRDS("./output/cod0_zinb_k3.rds")
cod0s_zinb_k3  <- readRDS("./output/cod0s_zinb_k3.rds")
cod1sg_zinb_k3 <- readRDS("./output/cod1sg_zinb_k3.rds")
cod2sg_zinb_k3 <- readRDS("./output/cod2sg_zinb_k3.rds")

loo(cod0_zinb_k3, cod0s_zinb_k3,
    cod1sg_zinb_k3, cod2sg_zinb_k3, moment_match = T, reloo = T)

looic <- c(cod0_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
           cod0s_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
           cod1sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
           cod2sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"])
df <- data.frame(model = c("cod0_zinb_k3", "cod0s_zinb_k3",
                           "cod1sg_zinb_k3", "cod2sg_zinb_k3"),
                 looic = looic)
df <- df[order(df$looic), ]
print(df)

# cod2sg_zinb_k3 is the selected model

## Predicted effects ---------------------------------------

## temp.anom predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod2sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod2sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod2sg_zinb_k3, effect = "temp.anom", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temp.anom
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temp.anom[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temp.anom[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temp.anom[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temp.anom[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(cod.data$temp.anom), amount = 0.1),
                          rep(NA, 100-length(unique(cod.data$temp.anom))))

fig.2a <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1, color = "red3") +
    labs(x = "Egg/larval temperature anomaly", y = "Fish / set") +
    scale_y_continuous(breaks=c(0,1,5,10,50,100,200)) +
    coord_trans(y = "pseudo_log") +
    theme_bw()+
    geom_rug(aes(x=rug.anom, y=NULL))
print(fig.2a)

ggsave("./figs/temp.anom_predicted_effect_cod2sg_zinb_k3.png", width = 3, height = 2)
## this is Fig. 2a in the draft


## Julian predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod2sg_zinb_k3, effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

g <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    coord_trans(y = "pseudo_log") +
    labs(x = "Day of year", y = "Cod abundance") +
    theme_bw()
print(g)
ggsave("./figs/julian_predicted_effect_cod2sg_zinb_k3.png", width = 5, height = 4)


## SSB predictions ##

## 95% CI
ce1s_1 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(cod2sg_zinb_k3, effect = "ssb", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$ssb
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$ssb[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$ssb[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$ssb[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$ssb[["lower__"]]

g <- ggplot(dat_ce) +
    aes(x = effect1__, y = estimate__) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
    geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
    geom_line(size = 1.5, color = "red3") +
    labs(x = "SSB", y = "Cod abundance") +
    coord_trans(y = "pseudo_log") +
    theme_bw()
print(g)
ggsave("./figs/SSB_predicted_effect_cod2sg3_zinb_k3.png", width = 5, height = 4)
