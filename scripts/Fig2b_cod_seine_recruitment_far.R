## parameterize cod cpue model with Fraction of Attibutable Risk as covariate
## this is Fig. 2b in the draft

library(dplyr)
library(plyr)
library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------


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

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

cod.data <- left_join(cod.data, obs)
cod.data$far_fac <- as.factor(if_else(cod.data$estimate__ >= 0.98, "high", "low"))

## brms: setup ---------------------------------------------

## Define model formula
cod_far_formula <-  bf(cod ~ s(julian, k = 4) + far_fac + s(ssb, k = 4) +
                        (1 | bay_fac/site_fac),
                      zi ~ s(julian, k = 4) + far_fac + s(ssb, k = 4) +
                        (1 | bay_fac/site_fac))


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

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
cod_far_zinb <- brm(cod_far_formula,
                    data = cod.data,
                    prior = priors_zinb,
                    family = zinb,
                    cores = 4, chains = 4, iter = 4000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 10))
cod_far_zinb  <- add_criterion(cod_far_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_far_zinb, file = "output/cod_far_zinb.rds")

cod_far_zinb <- readRDS("./output/cod_far_zinb.rds")
check_hmc_diagnostics(cod_far_zinb$fit)
neff_lowest(cod_far_zinb$fit)
rhat_highest(cod_far_zinb$fit)
summary(cod_far_zinb)
bayes_R2(cod_far_zinb)
plot(cod_far_zinb$criteria$loo, "k")
plot(conditional_smooths(cod_far_zinb), ask = FALSE)
y <- cod.data$cod
yrep_cod_far_zinb  <- fitted(cod_far_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_far_zinb[sample(nrow(yrep_cod_far_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_far_zinb")
pdf("./figs/trace_cod_far_zinb.pdf", width = 6, height = 4)
trace_plot(cod_far_zinb$fit)
dev.off()
            
## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_far_zinb, effect = "far_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$far_fac %>%
  select(far_fac, estimate__, lower__, upper__)

plot$far_fac <- reorder(plot$far_fac, desc(plot$far_fac))

fig.2b <- ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Fish / set") +
  xlab("FAR") +
  scale_x_discrete(labels=c(expression("<0.98"), expression("">=0.98))) +
  scale_y_continuous(breaks=c(1,5,10,50,100,150)) +
  coord_trans(y = "pseudo_log") + 
  theme_bw()

print(fig.2b)

ggsave("figs/prelim_FAR_recruit_plot.png", width=1.5, height=2, units='in')
## this is Fig. 2b
