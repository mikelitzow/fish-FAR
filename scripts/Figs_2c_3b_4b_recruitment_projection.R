## predict the response of cod and pollock recruitment to FAR values
## in the observational record and for CMIP projections
## currently used in Fig. 4b, and will become the new Figs. 2c, 3b

library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
library(reshape2)
source("./scripts/stan_utils.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())

## load data --------------------

# stock assessment model estimates
recr <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

## note that we're only using FAR as SSB wasn't important in original model

# 2017-2020 recruitment estimates are poorly supported by data; discard
recr$codR0.2020[recr$year >= 2017] <-  NA

recr <- recr %>%
  filter(year >= 1977) %>%
  select(year, codR0.2020) %>%
  mutate(sc.log.codR0=as.vector(scale(log(codR0.2020))))


# annual recruitment estimates from seines - load brms object
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")

seine <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))

print(seine)

seine$year_fac
seine.r <- data.frame(year=2006:2020,
                      seine.R=seine$year_fac$estimate__)

obs.dat <- left_join(recr, seine.r)

# log and scale seine estimates
# then predict 2017-2020 recruitment from seine values

obs.dat$sc.log.seine.R <- as.vector(scale(log(obs.dat$seine.R)))

mod <- lm(sc.log.codR0 ~ sc.log.seine.R, data=obs.dat, na.action = "na.omit")
summary(mod)

new.dat <- data.frame(sc.log.seine.R = obs.dat$sc.log.seine.R[obs.dat$year %in% 2017:2020])
obs.dat$sc.log.codR0[obs.dat$year %in% 2017:2020] <- predict(mod, newdata =  new.dat)

# finally, add FAR estimates for 1977:2020 - load brms object
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

obs.FAR <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))

obs.FAR <- obs.FAR$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__, se__)

names(obs.FAR)[2:3] <- c("FAR", "FAR.SE")

obs.dat <- left_join(obs.dat, obs.FAR)

## finally, load projected FAR values through 2046 - load brms object
mod_far_fixef <- readRDS("./output/mod_far_fixef.rds")

proj.far <- conditional_effects(mod_far_fixef, probs = c(0.025, 0.975))

print(proj.far)

proj.dat <- proj.far$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  filter(year >= 2006) %>% # limit to RCP8.5 projections (no historical observations)
  select(year, estimate__, se__)

names(proj.dat)[2:3] <- c("FAR", "FAR.SE")

## proj.dat is the dataframe for use as new.data when predicting 

## model recruitment as a function of FAR
R1 <- brm(sc.log.codR0 ~ s(FAR, k = 4),
          data = obs.dat,
          save_pars = save_pars(latent = TRUE),
          cores = 4, iter = 4000, chains = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(R1, file = "output/cod_R_FAR_obs.rds")
summary(R1)
names(R1$fit)


R1 <- readRDS("./output/cod_R_FAR_obs.rds")
bayes_R2(R1)
summary(R1)
plot(conditional_effects(R1), points = TRUE)
check_hmc_diagnostics(R1$fit)
neff_lowest(R1$fit)
rhat_highest(R1$fit)

# quick look at results
ce1s_1 <- conditional_effects(R1)
check <- ce1s_1$FAR %>%
  arrange(desc(FAR))
check

## try adding uncertainty in FAR
R2 <- brm(sc.log.codR0 ~ me(FAR, FAR.SE),
          data = obs.dat,
          save_pars = save_pars(latent = TRUE),
          cores = 4, iter = 4000, chains = 4,
          control = list(adapt_delta = 0.99, max_treedepth = 12))
saveRDS(R2, file = "output/cod_R_FAR_w_SE.rds")
summary(R2)
names(R2$fit)


R2 <- readRDS("./output/cod_R_FAR_w_SE.rds")
plot(conditional_effects(R2), points = TRUE)
check_hmc_diagnostics(R2$fit)
neff_lowest(R2$fit)
rhat_highest(R2$fit)

## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(R1, effect = "FAR", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(R1, effect = "FAR", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(R1, effect = "FAR", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$FAR
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$FAR[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$FAR[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$FAR[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$FAR[["lower__"]]

## and years for plot
cod.sub <- conditional_effects(R1, effect = "FAR")
cod.sub <- lapply(cod.sub, attributes)$FAR$points
cod.sub$year <- 1977:2020

# jitter x and y for plot
f <- 135 # set jitter factor
set.seed(22)
cod.sub$x.jitter <- jitter(cod.sub$FAR, factor=f)
cod.sub$y.jitter <- jitter(cod.sub$resp__, factor=f)

fig.2c <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  geom_text(data = cod.sub,
            aes(x = x.jitter, y = y.jitter, label = year), color = "grey40", size = 3) +
  labs(x = "Fraction of attributable risk", y = "Log recruitment anomaly") +
  theme_bw()
print(fig.2c)

ggsave("./figs/fig.2c.png", width = 4, height = 3)


# predict for CMIP projections!
post <- data.frame()
for(i in 2006:2046){
  # i <- 2006
  newdata <- data.frame(FAR=proj.dat$FAR[proj.dat$year==i],
                        FAR.SE=proj.dat$FAR.SE[proj.dat$year==i])
  
  xx <- data.frame(year = i, 
                   posterior = posterior_epred(R2, newdata = newdata, re_formula = NA, resp = "sc.log.codR0"))
  
  post <- rbind(post, xx)
  
}

post$decade <- ifelse(post$year %in% 2010:2019, "2010s",
                      ifelse(post$year %in% 2020:2029, "2020s",
                             ifelse(post$year %in% 2030:2039, "2030s", 
                                    ifelse(post$year <=2009, "2000s", "2040s"))))


# and add historical predictions
histor <- data.frame()
for(i in 1977:2019){
  # i <- 2006
  newdata <- data.frame(FAR=obs.FAR$FAR[obs.FAR$year==i],
                        FAR.SE=obs.FAR$FAR.SE[obs.FAR$year==i])
  
  xx <- data.frame(year = i, 
                   posterior = posterior_epred(R2, newdata = newdata, re_formula = NA, resp = "sc.log.codR0"))
  
  histor <- rbind(histor, xx)
  
}

histor$decade <- "Historical"

## combine and plot with CIs
pred.recr <- rbind(post, histor) %>%
  group_by(decade) %>%
  summarise(median=median(posterior),
            LCI=quantile(posterior, probs = 0.025),
            UCI=quantile(posterior, probs = 0.975)) %>%
  filter(decade %in% c("2020s", "2030s", "2040s", "Historical"))

pred.recr$order <- c(2,3,4,1)
pred.recr$decade <- reorder(pred.recr$decade, pred.recr$order)

cod.project.R <- ggplot(pred.recr, aes(decade, median)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2) +
  ylab("Log recruitment anomaly") +
  theme(axis.title.x = element_blank())

cod.project.R

ggsave("./figs/hist-projected_cod_R_with_FAR_uncertainty.png", width = 3, height = 3)

## now pollock! -----------------------
## load data --------------------

# stock assessment model estimates
recr <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

poll.obs.dat <- recr %>%
  filter(year %in% 1970:2019) %>%
  select(year, pollR0.2020, poll.SSB.2020) %>%
  mutate(sc.log.pollR0=as.vector(scale(log(pollR0.2020))))


obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

obs.FAR <- conditional_effects(obs_far_fixef, probs = c(0.025, 0.975))

obs.FAR <- obs.FAR$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  select(year, estimate__, se__)

names(obs.FAR)[2:3] <- c("FAR", "FAR.SE")

poll.obs.dat <- left_join(poll.obs.dat, obs.FAR)
names(poll.obs.dat)[3] <- "SSB"

## finally, load projected FAR values through 2046 - load brms object
mod_far_fixef <- readRDS("./output/mod_far_fixef.rds")

proj.far <- conditional_effects(mod_far_fixef, probs = c(0.025, 0.975))

print(proj.far)

proj.dat <- proj.far$year_fac %>%
  mutate(year = as.numeric(as.character(year_fac))) %>%
  filter(year >= 2006) %>% # limit to RCP8.5 projections (no historical observations)
  select(year, estimate__, se__)

names(proj.dat)[2:3] <- c("FAR", "FAR.SE")

## proj.dat is the dataframe for use as new.data when predicting 

## model recruitment as a function of FAR
poll.R1 <- brm(sc.log.pollR0 ~ s(FAR, k = 5),
          data = poll.obs.dat,
          seed = 1234,
          save_pars = save_pars(all = TRUE),
          cores = 4, iter = 4000, chains = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 12))
poll.R1 <- add_criterion(poll.R1, c("loo", "bayes_R2"),
                         moment_match = TRUE, reloo = FALSE,
                         cores = 4, k_threshold = 0.7)

saveRDS(poll.R1, file = "output/poll_R1_FAR_obs.rds")
summary(poll.R1)
bayes_R2(poll.R1)
names(poll.R1$fit)
plot(conditional_effects(poll.R1), points = TRUE)
check_hmc_diagnostics(poll.R1$fit)
neff_lowest(poll.R1$fit)
rhat_highest(poll.R1$fit)

## add SSB
# model selection
poll.R1s <- brm(sc.log.pollR0 ~ s(FAR, k = 5) + s(SSB, k = 5),
               data = poll.obs.dat,
               seed = 1234,
               save_pars = save_pars(all = TRUE),
               cores = 4, iter = 4000, chains = 4,
               control = list(adapt_delta = 0.999, max_treedepth = 12))
poll.R1s <- add_criterion(poll.R1s, c("loo", "bayes_R2"),
                          moment_match = TRUE, reloo = TRUE,
                          cores = 4, k_threshold = 0.7)

saveRDS(poll.R1s, file = "output/poll_R1s_FAR_obs.rds")
summary(poll.R1s)
names(poll.R1s$fit)
plot(conditional_effects(poll.R1s), points = TRUE)
check_hmc_diagnostics(poll.R1s$fit)
neff_lowest(poll.R1s$fit)
rhat_highest(poll.R1s$fit)


## model selection -----------------------------
poll.R1 <- readRDS("./output/poll_R1_FAR_obs.rds")
poll.R1s <- readRDS("./output/poll_R1s_FAR_obs.rds")

loo(poll.R1, poll.R1s)

## plot R1 for Fig. 3b ------------------------------------

## far predictions ##

## 95% CI
ce1s_1 <- conditional_effects(poll.R1, effect = "FAR", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(poll.R1, effect = "FAR", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(poll.R1, effect = "FAR", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$FAR
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$FAR[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$FAR[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$FAR[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$FAR[["lower__"]]
# dat_ce[["rug.anom"]] <- c(unique(poll.obs.dat$FAR),
#                           rep(NA, 100-length(unique(poll.obs.dat$FAR))))

fig.3b <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Fraction of attributable risk", y = "Log recruitment anomaly") +
  theme_bw()+
  geom_text(data = poll.obs.dat,
            aes(x = FAR, y = sc.log.pollR0, label = year), color = "grey40", size = 3) 
print(fig.3b)

## try adding uncertainty in FAR -----------------------------------------
poll.R2 <- brm(sc.log.pollR0 ~ 1 + me(FAR, FAR.SE) + I(me(FAR, FAR.SE)^2),
          data = poll.obs.dat,
          save_pars = save_pars(latent = TRUE),
          cores = 4, iter = 6000, chains = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(poll.R2, file = "output/poll_R_FAR_w_SE.rds")
summary(poll.R2)
names(poll.R2$fit)

poll.R3 <- brm(sc.log.pollR0 ~ 1 + me(FAR, FAR.SE) + I(me(FAR, FAR.SE)^2) + I(me(FAR, FAR.SE)^3),
               data = poll.obs.dat,
               save_pars = save_pars(latent = TRUE),
               cores = 4, iter = 6000, chains = 4,
               control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(poll.R3, file = "output/poll_R_FAR_w_SE_cubic.rds")
summary(poll.R3)
names(poll.R3$fit)


poll.R2 <- readRDS("./output/poll_R_FAR_w_SE.rds")

poll.R3 <- readRDS("./output/poll_R_FAR_w_SE_cubic.rds")

loo(poll.R2, poll.R3)

## Predict R2 manually
pred_far <- data.frame(FAR = seq(min(poll.obs.dat$FAR), max(poll.obs.dat$FAR), length.out = 100),
                       FAR.SE = 0.00001) ## set measurement error to zero
pred_full <- posterior_epred(poll.R2, newdata = pred_far)
pred <- data.frame(estimate = apply(pred_full, 2, mean),
                   upper = apply(pred_full, 2, quantile, probs = 0.975),
                   lower = apply(pred_full, 2, quantile, probs = 0.025))
pred_df <- cbind(pred_far, pred)
g <- ggplot(pred_df) +
    geom_ribbon(aes(x = FAR, ymin = lower, ymax = upper), fill = "grey90") +
    geom_line(aes(x = FAR, y = estimate), color = "red3") +
    geom_point(data = poll.obs.dat, aes(x = FAR, y = sc.log.pollR0), color = "grey25")
    ## geom_smooth(data = poll.obs.dat, aes(x = FAR, y = sc.log.pollR0), method = "lm",
    ##             formula = y ~ x + I(x^2), se = FALSE, color = "blue") +
    ## geom_segment(data = poll.obs.dat, aes(y = sc.log.pollR0, yend = sc.log.pollR0,
    ##                                  x = FAR - FAR.SE, xend = FAR + FAR.SE))
print(g)

check_hmc_diagnostics(poll.R2$fit)
neff_lowest(poll.R2$fit)
rhat_highest(poll.R2$fit)


## predict for CMIP projections!---------------------------
poll.post <- data.frame()
for(i in 2006:2046){
  # i <- 2006
  newdata <- data.frame(FAR=proj.dat$FAR[proj.dat$year==i],
                        FAR.SE=proj.dat$FAR.SE[proj.dat$year==i])
  
  xx <- data.frame(year = i, 
                   posterior = posterior_epred(poll.R3, newdata = newdata, re_formula = NA, resp = "sc.log.pollR0"))
  
  poll.post <- rbind(poll.post, xx)
  
}

poll.post$decade <- ifelse(poll.post$year %in% 2010:2019, "2010s",
                      ifelse(poll.post$year %in% 2020:2029, "2020s",
                             ifelse(poll.post$year %in% 2030:2039, "2030s", 
                                    ifelse(poll.post$year <=2009, "2000s", "2040s"))))


# and add historical (20th century) predictions
poll.histor <- data.frame()
for(i in 1977:2019){
  # i <- 2006
  newdata <- data.frame(FAR=obs.FAR$FAR[obs.FAR$year==i],
                        FAR.SE=obs.FAR$FAR.SE[obs.FAR$year==i])
  
  xx <- data.frame(year = i, 
                   posterior = posterior_epred(poll.R3, newdata = newdata, re_formula = NA, resp = "sc.log.codR0"))
  
  poll.histor <- rbind(poll.histor, xx)
  
}

poll.histor$decade <- "Historical"




## combine and plot with CIs
poll.pred.recr <- rbind(poll.post, poll.histor) %>%
  group_by(decade) %>%
  summarise(median=median(posterior),
            LCI=quantile(posterior, probs = 0.025),
            UCI=quantile(posterior, probs = 0.975)) %>%
  filter(decade %in% c("2020s", "2030s", "2040s", "Historical"))

poll.pred.recr$order <- c(2,3,4,1)
poll.pred.recr$decade <- reorder(poll.pred.recr$decade, poll.pred.recr$order)

ggplot(poll.pred.recr, aes(decade, median)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2) +
  ylab("Log recruitment anomaly") +
  theme(axis.title.x = element_blank())

ggsave("./figs/hist-projected_poll_R_with_FAR_uncertainty.png", width = 3, height = 3)

# combine projections into a single plot

pred.recr$species <- "cod"
poll.pred.recr$species <- "pollock"

all.plot <- rbind(pred.recr, poll.pred.recr)

cod.poll.proj.R <- ggplot(all.plot, aes(decade, median, color=species)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.2, position=position_dodge(width=0.5)) +
  ylab("Log recruitment anomaly") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.8)) +
  scale_color_manual(values=cb[c(2,4)])

cod.poll.proj.R

ggsave("./figs/hist-projected_poll_cod_ R.png", width = 3, height = 3)
## this is Fig. 4b in the draft

## calculate change in medians ------------------------------------
summ.dat <- all.plot %>%
  filter(decade %in% c("2020s", "Historical")) %>%
  select(species, decade, median) %>%
  pivot_wider(values_from = median, names_from = species)

# go back to original data frames to back-calculate to original units!     
ggplot(obs.dat, aes(sc.log.codR0, codR0.2020)) +
  geom_point()

cod.mod <- lm(log(codR0.2020) ~ sc.log.codR0, data=obs.dat, na.action = "na.exclude")
summary(cod.mod)

check <- exp(predict(cod.mod))

plot(check, obs.dat$codR0.2020) # right!
new.dat <- data.frame(sc.log.codR0 = summ.dat$cod)
summ.dat$raw.cod <- exp(predict(cod.mod, newdata = new.dat))

# and pollock
poll.mod <- lm(log(pollR0.2020) ~ sc.log.pollR0, data=poll.obs.dat, na.action = "na.exclude")
summary(poll.mod)

check <- exp(predict(poll.mod))

plot(check, poll.obs.dat$pollR0.2020) # right!
new.dat <- data.frame(sc.log.pollR0 = summ.dat$pollock)
summ.dat$raw.poll <- exp(predict(poll.mod, newdata = new.dat))

cod.change <- (summ.dat$raw.cod[2]-summ.dat$raw.cod[1])/summ.dat$raw.cod[2]
poll.change <- (summ.dat$raw.poll[2]-summ.dat$raw.poll[1])/summ.dat$raw.poll[2]
