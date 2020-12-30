## predict the response of cod and pollock recruitment to FAR values
## in the observational record and for CMIP projections
## currently used in Fig. 4b, and will become the new Figs. 2c, 3c

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

# 2018-2020 recruitment estimates are poorly supported by data; discard
recr$codR0.2020[recr$year >= 2018] <-  NA

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
# then predict 2018-2020 recruitment from seine values

obs.dat$sc.log.seine.R <- as.vector(scale(log(obs.dat$seine.R)))

mod <- lm(sc.log.codR0 ~ sc.log.seine.R, data=obs.dat, na.action = "na.omit")
summary(mod)

new.dat <- data.frame(sc.log.seine.R = obs.dat$sc.log.seine.R[obs.dat$year %in% 2018:2020])
obs.dat$sc.log.codR0[obs.dat$year %in% 2018:2020] <- predict(mod, newdata =  new.dat)

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
plot(conditional_effects(R1), points = TRUE)

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

## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(R2, effect = "FAR", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(R2, effect = "FAR", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(R2, effect = "FAR", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$FAR
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$FAR[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$FAR[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$FAR[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$FAR[["lower__"]]

## and years for plot
cod.sub <- conditional_effects(R2, effect = "FAR")
cod.sub <- lapply(cod.sub, attributes)$FAR$points
cod.sub$year <- 1977:2020

cod_R_FAR_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  geom_text(data = cod.sub,
            aes(x = FAR, y = resp__, label = year), color = "grey40", size = 3) +
  labs(x = "Fraction of attributable risk", y = "Recruitment anomaly") +
  theme_bw()
print(cod_R_FAR_plot)

ggsave("./figs/predicted_effect_cod_R_FAR_with_uncertainty.png", width = 4, height = 3)
## this is Fig. 2c in the draft

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

##-------------------------------
## now pollock!
## load data --------------------

# stock assessment model estimates
recr <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

# need to estimate pre-1977 recruitment,
# which was provided in the 2019 report but not in 2020
mod <- lm(pollR0.2020 ~ pollR0.2019, data=recr)
summary(mod)
pred.R0 <- predict(mod, newdata = recr[recr$year < 1977,])
recr$pollR0.2020[recr$year < 1977] <- pred.R0

# 2020 recruitment estimate is poorly supported by data; discard
# note that this differs from cod model! better support for recent year classes in pollock
# from the acoustic survey
recr$pollR0[recr$year == 2020] <-  NA

recr <- recr %>%
  filter(year >= 1970) %>%
  select(year, pollR0.2020) %>%
  mutate(sc.log.pollR0=as.vector(scale(log(pollR0.2020))))


# add 2020 as NA
xtra <- data.frame(year=2020,
                   pollR0.2020=NA,
                   sc.log.pollR0=NA)

recr <- rbind(recr, xtra)

# annual recruitment estimates from seines - load brms object
recr_2_zinb <- readRDS("./output/poll_recr_2_zinb.rds")

seine <- conditional_effects(recr_2_zinb, effect = "year_fac", re_formula = NA,
                             probs = c(0.025, 0.975))

print(seine)

seine$year_fac
seine.r <- data.frame(year=2006:2020,
                      seine.R=seine$year_fac$estimate__)

obs.dat <- left_join(recr, seine.r)

# log and scale seine estimates
# then predict 2020 recruitment from seine values

obs.dat$sc.log.seine.R <- as.vector(scale(log(obs.dat$seine.R)))

mod <- lm(sc.log.pollR0 ~ sc.log.seine.R, data=obs.dat, na.action = "na.omit")
summary(mod)

new.dat <- data.frame(sc.log.seine.R = obs.dat$sc.log.seine.R[obs.dat$year == 2020])
obs.dat$sc.log.pollR0[obs.dat$year == 2020] <- predict(mod, newdata =  new.dat)

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
poll.R1 <- brm(sc.log.pollR0 ~ s(FAR, k = 4),
          data = obs.dat,
          save_pars = save_pars(latent = TRUE),
          cores = 4, iter = 4000, chains = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 12))
saveRDS(poll.R1, file = "output/poll_R_FAR_obs.rds")
summary(poll.R1)
names(poll.R1$fit)


poll.R1 <- readRDS("./output/poll_R_FAR_obs.rds")
plot(conditional_effects(poll.R1), points = TRUE)

## try adding uncertainty in FAR
poll.R2 <- brm(sc.log.pollR0 ~ 1 + me(FAR, FAR.SE) + I(me(FAR, FAR.SE)^2),
          data = obs.dat,
          save_pars = save_pars(latent = TRUE),
          cores = 4, iter = 6000, chains = 4,
          control = list(adapt_delta = 0.999, max_treedepth = 16))

saveRDS(poll.R2, file = "output/poll_R_FAR_w_SE.rds")
summary(poll.R2)
names(poll.R2$fit)


poll.R2 <- readRDS("./output/poll_R_FAR_w_SE.rds")
plot(conditional_effects(poll.R2), points = TRUE) # this isn't working!

check_hmc_diagnostics(poll.R2$fit)
neff_lowest(R2$fit)
rhat_highest(R2$fit)


## predict for CMIP projections!---------------------------
poll.post <- data.frame()
for(i in 2006:2046){
  # i <- 2006
  newdata <- data.frame(FAR=proj.dat$FAR[proj.dat$year==i],
                        FAR.SE=proj.dat$FAR.SE[proj.dat$year==i])
  
  xx <- data.frame(year = i, 
                   posterior = posterior_epred(poll.R2, newdata = newdata, re_formula = NA, resp = "sc.log.pollR0"))
  
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
                   posterior = posterior_epred(poll.R2, newdata = newdata, re_formula = NA, resp = "sc.log.codR0"))
  
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
  scale_color_manual(values=cb[3:4])

cod.poll.proj.R

ggsave("./figs/hist-projected_poll_cod_ R.png", width = 3, height = 3)
## this is Fig. 4b in the draft