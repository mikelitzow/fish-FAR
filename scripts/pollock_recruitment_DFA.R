## estimate age-0 pollock abundance from multiple data sources

library(tidyverse)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./scripts/stan_utils.R")
library(MARSS)
theme_set(theme_bw())

## Fit brms model to beach seine data --------------------------------------------
# read in data
poll.data <- read.csv("data/cpue.data.csv")
poll.data$bay_fac <- as.factor(poll.data$bay)
poll.data$year_fac <- as.factor(poll.data$year)
poll.data$site_fac <- as.factor(poll.data$site)
poll.data$bay_site_fac <- as.factor(paste0(poll.data$bay, "_", poll.data$site))
poll.data$present <- ifelse(poll.data$pollock > 0, 1, 0)
poll.data$date <- as.Date(poll.data$julian,
                          origin = paste0(poll.data$year, "-01-01"))
# restrict to long-term sites and AK Peninsula bays with high positive catches
levels(poll.data$bay_fac)
keep <- c("Agripina", "Anton Larson Bay", "Balboa", "Cook Bay", "Mitrofania", "Port Wrangell") 
poll.data <- poll.data %>%
  filter(bay_fac %in% keep)

## Check distributions
plot(poll.data$pollock) # even more zeros than cod, as expected
hist(poll.data$pollock, breaks = 100) ## lots of zeros
tab <- table(poll.data$pollock)
plot(tab)
summary(stats::glm(pollock ~ 1, data = poll.data, family = poisson))
summary(MASS::glm.nb(pollock ~ 1, data = poll.data))

## Percent zeros by bay
plyr::ddply(poll.data, .(bay), summarize,
            zeros = sum(pollock == 0),
            not_zeros = sum(pollock > 0),
            perc_zero = (zeros / (zeros + not_zeros)) * 100)


g <- ggplot(poll.data) +
  aes(x = date, y = pollock, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

## brms: setup ---------------------------------------------

## Define model formula

recr_2_formula <-  bf(pollock ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac),
                      zi ~ s(julian, k = 3) + (1 | bay_fac/site_fac) + (year_fac))
# (this is the best model from the full seine data set [all bays])

## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Show default priors
get_prior(recr_2_formula, poll.data, family = zinb)


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

poll_recr_2_zinb_reduced_bays <- brm(recr_2_formula,
                        data = poll.data,
                        prior = priors_zinb,
                        family = zinb,
                        cores = 4, chains = 4, iter = 6000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.99, max_treedepth = 10))
poll_recr_2_zinb_reduced_bays  <- add_criterion(poll_recr_2_zinb_reduced_bays, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(poll_recr_2_zinb_reduced_bays, file = "output/poll_recr_2_zinb_reduced_bays.rds")

poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")
check_hmc_diagnostics(poll_recr_2_zinb_reduced_bays$fit)
neff_lowest(poll_recr_2_zinb_reduced_bays$fit)
rhat_highest(poll_recr_2_zinb_reduced_bays$fit)
summary(poll_recr_2_zinb_reduced_bays)
bayes_R2(poll_recr_2_zinb_reduced_bays)
plot(poll_recr_2_zinb_reduced_bays$criteria$loo, "k")
plot(conditional_smooths(poll_recr_2_zinb_reduced_bays), ask = FALSE)
y <- poll.data$pollock
yrep_poll_recr_2_zinb_reduced_bays  <- fitted(poll_recr_2_zinb_reduced_bays, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_poll_recr_2_zinb_reduced_bays[sample(nrow(yrep_poll_recr_2_zinb_reduced_bays), 25), ]) +
  xlim(0, 500) +
  ggtitle("poll_recr_2_zinb_reduced_bays")
pdf("./figs/trace_poll_recr_2_zinb_reduced_bays.pdf", width = 6, height = 4)
trace_plot(poll_recr_2_zinb_reduced_bays$fit)
dev.off()

# extract (and plot) annual estimates
poll_recr_2_zinb_reduced_bays <- readRDS("./output/poll_recr_2_zinb_reduced_bays.rds")

ce1s_1 <- conditional_effects(poll_recr_2_zinb_reduced_bays, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(ce1s_1)
ggsave("./figs/annual_poll_recruitment_estimates_poll_recr_2_zinb_reduced_bays.png", width = 7, height = 4)

# load eco-foci larval / age-0 abundance data
foci.larv <- read.csv("data/ECO-FOCI_larval_pollock.csv")
head(foci.larv)

foci.juv <- read.csv("data/ECO-FOCI_age_0_pollock.csv")
head(foci.juv)

# load mace age-1 index
mace <- read.csv("./data/MACE_age_1_pollock.csv")
head(mace)
mace$year.0 <- mace$year-1

# combine the four data sets
seine.dat <- ce1s_1$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__)

names(seine.dat)[2] <- "seine.est"
names(foci.larv)[2:3] <- c("larv.est", "year") 
names(foci.juv)[1:2] <- c("year", "juv.est")

mace <- mace %>%
  select(-year)
names(mace) <- c("mace.est.lag1", "year")

dat <- data.frame(year = 1981:2020)

dat <- left_join(dat, foci.larv)
dat <- left_join(dat, foci.juv)
dat <- left_join(dat, mace)
dat <- left_join(dat, seine.dat)

head(dat)

# clean up!
dat <- dat %>%
  select(year, larv.est, juv.est, mace.est.lag1, seine.est)

# now log-transform and scale!
scaled.dat <- dat
for(j in 2:5){
  
  scaled.dat[,j] <- as.vector(scale(log(dat[,j]+0.01)))
  # adding a constant b/c mace has 1 year = 0!
  
}

# check correlations
cor(scaled.dat[,2:5], use="p") # pretty strong!


## fit a DFA model ---------------------------------------------

# set up data
dfa.dat <- as.matrix(t(scaled.dat[,2:5]))
colnames(dfa.dat) <- scaled.dat$year

## compare with reduced model (no MACE data, which is alrady in assessment!)
reduced.dat <- dfa.dat[rownames(dfa.dat) != "mace.est.lag1",]

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # allowing up to 1 trends
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(reduced.dat[,colnames(reduced.dat) %in% 1987:2020], model=dfa.model,
                 form="dfa", z.score=TRUE, control=cntl.list)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

## best model is equalvarcov - but that returns loadings of 0!
# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.list = list(A="zero", m=1, R="diagonal and equal")
dfa.mod = MARSS(reduced.dat[,colnames(reduced.dat) %in% 1987:2020], model=model.list, z.score=TRUE, form="dfa")

# get CI and plot loadings...
modCI <- MARSSparamCIs(dfa.mod)
modCI

loadings <- data.frame(names = c("Larval", "Age-0 trawl", "Age-0 seine"),
                       loading = modCI$par$Z,
                       upCI = modCI$par.upCI$Z,
                       lowCI = modCI$par.lowCI$Z)

ggplot(loadings, aes(names, loading)) +
  geom_bar(stat="identity", fill="light grey") +
  geom_errorbar(aes(ymin=lowCI, ymax=upCI), width=0.2) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.5)) 

ggsave("./figs/reduced_pollock_DFA_loadings.png", width=2, height=2, units = 'in')


# plot trend
trend <- data.frame(year = 1987:2020,
                    trend = as.vector(dfa.mod$states),
                    ymin = as.vector(dfa.mod$states-1.96*dfa.mod$states.se),
                    ymax = as.vector(dfa.mod$states+1.96*dfa.mod$states.se))

ggplot(trend, aes(year, trend)) +
  geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="grey90") +
  geom_line(color="red") 

write.csv(trend, "./output/poll_dfa_trend.csv")

ggsave("./figs/reduced_pollock_DFA_trend.png", width=3, height=2, units = 'in')

## fit brms model - FAR as factor ------------------------------------

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

trend <- left_join(trend, obs)

trend$far_fac <- as.factor(if_else(trend$estimate__ >= 0.95, "high", "low"))


## brms: setup ---------------------------------------------

## Define model formula
poll_dfa_far_formula <-  bf(trend ~  far_fac)

poll_dfa_far <- brm(poll_dfa_far_formula,
                    data = trend,
                    cores = 4, chains = 4, iter = 4000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(poll_dfa_far, file = "output/poll_dfa_far.rds")

poll_dfa_far <- readRDS("./output/poll_dfa_far.rds")
check_hmc_diagnostics(poll_dfa_far$fit)
neff_lowest(poll_dfa_far$fit)
rhat_highest(poll_dfa_far$fit)
summary(poll_dfa_far)
bayes_R2(poll_dfa_far)

## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(poll_dfa_far, effect = "far_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$far_fac %>%
  select(far_fac, estimate__, lower__, upper__)

plot$far_fac <- reorder(plot$far_fac, desc(plot$far_fac))

ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Recruitment anomaly") +
  xlab("FAR") +
  scale_x_discrete(labels=c(expression(""<=0.9), expression("">=0.98))) +
  theme_bw() 


ggsave("figs/prelim_FAR_poll_dfa_plot.png", width=1.5, height=2, units='in')

## Fit brms model - FAR as continuous variable -------------------
# add 2020 SSB values
ssb <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

ssb <- ssb %>%
  select(year, poll.SSB.2020)

names(ssb)[2] <- "ssb"

trend <- left_join(trend, ssb)

trend <- trend %>%
  select(estimate__, trend, ssb, year)

names(trend)[1] <- "far"

## plot data
g <- ggplot(trend) +
  aes(x = far, y = trend) +
  geom_point()
print(g)

g <- ggplot(trend) +
  aes(x = ssb, y = trend) +
  geom_point()
print(g)


## brms: setup ---------------------------------------------

## Define model formulas
## Limiting knots to 3 to prevent overfitting

dfa1_far_formula <-  bf(trend ~ s(ssb, k = 3) + s(far, k = 3))

dfa2_far_formula <-  bf(trend ~ s(far, k = 3))


## fit --------------------------------------
dfa1_far_brm <- brm(dfa1_far_formula,
                data = trend,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99, max_treedepth = 10))
dfa1_far_brm  <- add_criterion(dfa1_far_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(dfa1_far_brm, file = "output/dfa1_far_brm.rds")

dfa1_far_brm <- readRDS("./output/dfa1_far_brm.rds")
check_hmc_diagnostics(dfa1_far_brm$fit)
neff_lowest(dfa1_far_brm$fit)
rhat_highest(dfa1_far_brm$fit)
summary(dfa1_far_brm)
bayes_R2(dfa1_far_brm)
plot(dfa1_far_brm$criteria$loo, "k")
plot(conditional_smooths(dfa1_far_brm), ask = FALSE)
y <- trend$trend
yrep_dfa1_far_brm  <- fitted(dfa1_far_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_dfa1_far_brm[sample(nrow(yrep_dfa1_far_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("dfa1_far_brm")
pdf("./figs/trace_dfa1_far_brm.pdf", width = 6, height = 4)
trace_plot(dfa1_far_brm$fit)
dev.off()


dfa_2_far_brm <- brm(dfa2_far_formula,
                data = trend,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99, max_treedepth = 10))
dfa_2_far_brm  <- add_criterion(dfa_2_far_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(dfa_2_far_brm, file = "output/dfa_2_far_brm.rds")

dfa_2_far_brm <- readRDS("./output/dfa_2_far_brm.rds")
check_hmc_diagnostics(dfa_2_far_brm$fit)
neff_lowest(dfa_2_far_brm$fit)
rhat_highest(dfa_2_far_brm$fit)
summary(dfa_2_far_brm)
bayes_R2(dfa_2_far_brm)
plot(dfa_2_far_brm$criteria$loo, "k")
plot(conditional_smooths(dfa_2_far_brm), ask = FALSE)
y <- trend$trend
yrep_dfa_2_far_brm  <- fitted(dfa_2_far_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_dfa_2_far_brm[sample(nrow(yrep_dfa_2_far_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("dfa_2_far_brm")
pdf("./figs/trace_dfa_2_far_brm.pdf", width = 6, height = 4)
trace_plot(dfa_2_far_brm$fit)
dev.off()


## Model selection -----------------------------------------
dfa1_far_brm  <- readRDS("./output/dfa1_far_brm.rds")
dfa2_far_brm  <- readRDS("./output/dfa_2_far_brm.rds")

loo(dfa1_far_brm, dfa2_far_brm)

## plot predicted values ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(dfa1_far_brm, effect = "far", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(dfa1_far_brm, effect = "far", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(dfa1_far_brm, effect = "far", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$far
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$far[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$far[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$far[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$far[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(trend$far), amount = 0.01),
                         rep(NA, 100-length(unique(trend$far))))


fig.3a <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fraction of attributable risk", y = "DFA trend") +
  theme_bw()+
 geom_rug(aes(x=rug.anom, y=NULL))
print(fig.3a)

ggsave("./figs/continuous_far_predicted_effect_dfa1_far_brm.png", width = 3, height = 2)

ce1s_1$far
## aside - how well does the DFA trend predict model-estimated recruitment?---------------

## add modeled R
r <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

r <- r %>%
  select(year, pollR0.2020)

names(r)[2] <- "R0"

trend <- left_join(trend, r)


## fit --------------------------------------
dfa1_far_formula <-  bf(trend ~ s(ssb, k = 3) + s(far, k = 3))

dfa_model_R_brm <- brm(bf(log(R0) ~ trend),
                    data = trend,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(dfa_model_R_brm, file = "output/dfa_model_R_brm.rds")

check_hmc_diagnostics(dfa_model_R_brm$fit)
neff_lowest(dfa_model_R_brm$fit)
rhat_highest(dfa_model_R_brm$fit)
summary(dfa_model_R_brm)
bayes_R2(dfa_model_R_brm)

y <- trend$trend
yrep_dfa_model_R_brm  <- fitted(dfa_model_R_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_dfa_model_R_brm[sample(nrow(yrep_dfa_model_R_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("dfa_model_R_brm")
pdf("./figs/trace_dfa_model_R_brm.pdf", width = 6, height = 4)
trace_plot(dfa_model_R_brm$fit)
dev.off()

## plot predicted values ---------------------------------------

dfa_model_R_brm  <- readRDS("./output/dfa_model_R_brm.rds")

## 95% CI
ce1s_1 <- conditional_effects(dfa_model_R_brm, effect = "trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(dfa_model_R_brm, effect = "trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(dfa_model_R_brm, effect = "trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$trend[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(trend$trend), amount = 0.01),
                          rep(NA, 100-length(unique(trend$trend))))

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "ln(model recruitment)") +
  theme_bw()+
  geom_rug(aes(x=rug.anom, y=NULL))
print(g)

ggsave("./figs/predicted_model_R_dfa_brm.png", width = 3, height = 2)
## Fit brms model - temp anomalies -----------------------------

# load GODAS anomalies
temp <- read.csv("./data/pollock godas anomalies.csv")

# check for correlation
cor(temp$egg[temp$year >=1987], temp$larval[temp$year >=1987])

# r = 0.53 - likely too high to treat separately!
temp <- temp %>%
  select(year, mean.anom)

trend <- left_join(trend, temp)

# add 2020 SSB values
ssb <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

ssb <- ssb %>%
  select(year, poll.SSB.2020)

names(ssb)[2] <- "ssb"

trend <- left_join(trend, ssb)

## plot data
g <- ggplot(trend) +
  aes(x = mean.anom, y = trend) +
  geom_point()
print(g)

g <- ggplot(trend) +
  aes(x = ssb, y = trend) +
  geom_point()
print(g)

## brms: setup ---------------------------------------------

## Define model formulas
## Limiting knots to 3 to prevent overfitting

dfa1_formula <-  bf(trend ~ s(ssb, k = 3) + s(mean.anom, k = 3))

dfa2_formula <-  bf(trend ~ s(mean.anom, k = 3))


## fit --------------------------------------
dfa1_brm <- brm(dfa1_formula,
                    data = trend,
                    cores = 4, chains = 4, iter = 3000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))
dfa1_brm  <- add_criterion(dfa1_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(dfa1_brm, file = "output/dfa1_brm.rds")

dfa1_brm <- readRDS("./output/dfa1_brm.rds")
check_hmc_diagnostics(dfa1_brm$fit)
neff_lowest(dfa1_brm$fit)
rhat_highest(dfa1_brm$fit)
summary(dfa1_brm)
bayes_R2(dfa1_brm)
plot(dfa1_brm$criteria$loo, "k")
plot(conditional_smooths(dfa1_brm), ask = FALSE)
y <- trend$trend
yrep_dfa1_brm  <- fitted(dfa1_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_dfa1_brm[sample(nrow(yrep_dfa1_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("dfa1_brm")
pdf("./figs/trace_dfa1_brm.pdf", width = 6, height = 4)
trace_plot(dfa1_brm$fit)
dev.off()


dfa2_brm <- brm(dfa2_formula,
                data = trend,
                cores = 4, chains = 4, iter = 3000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.99, max_treedepth = 10))
dfa2_brm  <- add_criterion(dfa2_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(dfa2_brm, file = "output/dfa2_brm.rds")

dfa2_brm <- readRDS("./output/dfa2_brm.rds")
check_hmc_diagnostics(dfa2_brm$fit)
neff_lowest(dfa2_brm$fit)
rhat_highest(dfa2_brm$fit)
summary(dfa2_brm)
bayes_R2(dfa2_brm)
plot(dfa2_brm$criteria$loo, "k")
plot(conditional_smooths(dfa2_brm), ask = FALSE)
y <- trend$trend
yrep_dfa2_brm  <- fitted(dfa2_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_dfa2_brm[sample(nrow(yrep_dfa2_brm), 25), ]) +
  xlim(0, 500) +
  ggtitle("dfa2_brm")
pdf("./figs/trace_dfa2_brm.pdf", width = 6, height = 4)
trace_plot(dfa2_brm$fit)
dev.off()


## Model selection -----------------------------------------
dfa1_brm  <- readRDS("./output/dfa1_brm.rds")
dfa2_brm  <- readRDS("./output/dfa2_brm.rds")

loo(dfa1_brm, dfa2_brm)

## DFA as a predictor of model R - brms -------------------------
# load assessment time series
recr <- read.csv("data/cod_pollock_assessment_2020_SAFEs.csv", row.names = 1)

# limiting to 2019 and earlier 
plot <- data.frame(year=1987:2019,
                   ln_assessment_model_R=log(recr$pollR0.2020[row.names(recr) %in% 1987:2019]),
                   dfa_trend=trend$trend[trend$year %in% 1987:2019])

cor(plot) # r = 0.64

ggplot(plot, aes(dfa_trend, ln_assessment_model_R)) +
  geom_text(aes(label=year)) +
  ggtitle("1987-2019") +
  theme_bw()

## fit a brms model ------------------------------------------

# first, clean up data
dat <- plot

names(dat)[2] <- "model"

## brms: setup ---------------------------------------------

## Define model formulas
## Limiting knots to 3 to prevent overfitting

pollR1_formula <-  bf(model ~ s(dfa_trend, k = 3))


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


## plot predicted values pollR1_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(pollR1_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR1_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR1_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa_trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa_trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa_trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa_trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa_trend[["lower__"]]

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "ln(model recruitment)") +
  geom_text(data=dat, aes(dfa_trend, model, label = year), size=1.5) +
  theme_bw()

print(g)

ggsave("./figs/pollR1_brm_dfa.png", width=3, height=2.5, units = 'in')

## second round of model-fitting adding seine:FAR interaction ------------

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

ce1s_1 <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))
obs <- ce1s_1$year_fac %>%
  select(year_fac, estimate__)

obs$year <- as.numeric(as.character(obs$year_fac))

dat <- left_join(dat, obs)
names(dat)[5] <- "far"

## Define model formulas

pollR2_formula <-  bf(model ~ s(dfa_trend, k = 3) + s(far, k=3))


## fit --------------------------------------
pollR2_brm <- brm(pollR2_formula,
                 data = dat,
                 cores = 4, chains = 4, iter = 3000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))
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

## plot predicted values pollR2_brm ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(pollR2_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR2_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR2_brm, effect = "dfa_trend", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$dfa_trend
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$dfa_trend[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$dfa_trend[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$dfa_trend[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$dfa_trend[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$dfa_trend), amount = 0.01),
                          rep(NA, 100-length(unique(dat$dfa_trend))))

g <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "DFA trend", y = "ln(model recruitment)") +
  theme_bw()+
  geom_rug(aes(x=rug.anom, y=NULL))
print(g)

ggsave("./figs/pollR2_brm_dfa.png", width=3, height=2, units = 'in')


## and far
## plot predicted values ---------------------------------------
## 95% CI
ce1s_1 <- conditional_effects(pollR2_brm, effect = "far", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(pollR2_brm, effect = "far", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(pollR2_brm, effect = "far", re_formula = NA,
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

ggsave("./figs/pollR2_brm_far.png", width=3, height=2, units = 'in')


## Bayes DFA - not very successful! -----------------------------
dfa.dat <- as.matrix(t(scaled.dat[,2:5]))
colnames(dfa.dat) <- scaled.dat$year

dfa.mod <- fit_dfa(y = dfa.dat[,colnames(dfa.dat) %in% 2000:2020], iter = 4000) 

str(dfa.mod)

print(plot_loadings(rotated, names = names))

plot_fitted(dfa.mod)

rotate <- rotate_trends(dfa.mod)
plot_trends(rotate, years = 2000:2020)
plot_loadings(rotate)
rownames(dfa.dat)
