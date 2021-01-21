library(tidyverse)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_bw())


## far panel ---------

# load FAR estimates
obs_far_fixef <- readRDS("./output/obs_far_fixef.rds")

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

far <- ggplot(filter(pred.obs, year >= 1969)) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  scale_x_continuous(breaks=seq(1960, 2020, 10)) 

print(far)



# and combine!
d1 <- d1 %>%
  select(year, mean.anom)
names(d1)[2] <- "Pollock"

d2 <- d2 %>%
  select(year, mean.anom)
names(d2)[2] <- "Cod"

dat <- left_join(d1, d2)

dat <- left_join(obs, dat)

dat <- dat %>%
  pivot_longer(cols = c(-year, -UCI, -LCI))

dat$panel = if_else(dat$name == "FAR", "Fraction of attributable risk", "Temperature anomaly (SD)")

ggplot(filter(dat, panel=="Temperature anomaly (SD)"), aes(year, value, color=name)) +
  geom_line() 


## model SSB and R estimates --------------------------------------------------------------------

mod.dat <- read.csv("./data/cod_pollock_assessment_2020_SAFEs.csv")

# clean up 
poll.r <- mod.dat %>%
  select(year, pollR0.2020, pollR0.2020.LCI, pollR0.2020.UCI) 
names(poll.r)[2:4] <- c("est", "LCI", "UCI")  
poll.r$variable <- "recruitment"

poll.s <- mod.dat %>%
  select(year, poll.SSB.2020, poll.SSB.2020.LCI, poll.SSB.2020.UCI) 
names(poll.s)[2:4] <- c("est", "LCI", "UCI")  
poll.s$variable <- "ssb"

poll <- rbind(poll.r, poll.s)
poll$sp <- "pollock"

cod.r <- mod.dat %>%
  select(year, codR0.2020, codR0.2020.SD) %>%
  mutate(LCI=codR0.2020-2*codR0.2020.SD, UCI=codR0.2020+2*codR0.2020.SD) %>%
  select(-codR0.2020.SD)

names(cod.r)[2] <- "est"  
cod.r$variable <- "recruitment"

cod.s <- mod.dat %>%
  select(year, codSSB.2020, codSSB.2020.SD) %>%
  mutate(LCI=codSSB.2020-2*codSSB.2020.SD, UCI=codSSB.2020+2*codSSB.2020.SD) %>%
  select(-codSSB.2020.SD)

names(cod.s)[2] <- "est"  
cod.s$variable <- "ssb"

cod <- rbind(cod.r, cod.s)
cod$sp <- "cod"

plot.dat <- rbind(cod, poll)

# change into correct units!
plot.true <- plot.dat
plot.true$est <- case_when(
  
  plot.dat$sp == "cod" & plot.dat$variable == "recruitment" ~ plot.dat$est * 1e+09,
  plot.dat$sp == "pollock" & plot.dat$variable == "recruitment" ~ plot.dat$est * 1e+06,
  
  plot.dat$sp == "cod" & plot.dat$variable == "ssb" ~ plot.dat$est,
  plot.dat$sp == "pollock" & plot.dat$variable == "ssb" ~ plot.dat$est * 1e+03,
  
)

plot.true$LCI <- case_when(
  
  plot.dat$sp == "cod" & plot.dat$variable == "recruitment" ~ plot.dat$LCI * 1e+09,
  plot.dat$sp == "pollock" & plot.dat$variable == "recruitment" ~ plot.dat$LCI * 1e+06,
  
  plot.dat$sp == "cod" & plot.dat$variable == "ssb" ~ plot.dat$LCI,
  plot.dat$sp == "pollock" & plot.dat$variable == "ssb" ~ plot.dat$LCI * 1e+03,
  
)

plot.true$UCI <- case_when(
  
  plot.dat$sp == "cod" & plot.dat$variable == "recruitment" ~ plot.dat$UCI * 1e+09,
  plot.dat$sp == "pollock" & plot.dat$variable == "recruitment" ~ plot.dat$UCI * 1e+06,
  
  plot.dat$sp == "cod" & plot.dat$variable == "ssb" ~ plot.dat$UCI,
  plot.dat$sp == "pollock" & plot.dat$variable == "ssb" ~ plot.dat$UCI * 1e+03,
  
)


plot.s <- plot.true %>%
  filter(variable == "ssb")

ssb.fig <- ggplot(plot.s, aes(year, est/1000, color=sp, fill=sp)) +
  geom_line() +
  coord_trans(y = "pseudo_log") +
  geom_ribbon(aes(ymin=UCI/1000, ymax=LCI/1000), alpha = 0.2, lty=0) +
  scale_color_manual(values=cb[3:4]) +
  scale_fill_manual(values=cb[3:4]) +
  scale_y_continuous(breaks=c(50, 100, 200, 400, 600, 800)) +
  ylab("Thousands of tons") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.4, 0.2))


ssb.fig


plot.r <- plot.true %>%
  filter(variable == "recruitment")

# and drop 2017-2020 cod!
drop <- plot.r$sp == "cod" & plot.r$variable == "recruitment" & plot.r$year %in% 2017:2020
plot.r <- plot.r[!drop,]

plot.r$LCI[which(plot.r$LCI < 0)] <- 1e+07

r.fig <- ggplot(plot.r, aes(year, est/1e+06, color=sp, fill=sp)) +
  geom_line() +
  facet_wrap(~sp, scales="free_y", nrow=2) +
  coord_trans(y = "pseudo_log") +
  geom_ribbon(aes(ymin=UCI/1e+06, ymax=LCI/1e+06), alpha = 0.2, lty=0) +
  scale_color_manual(values=cb[3:4]) +
  scale_fill_manual(values=cb[3:4]) +
  scale_y_continuous(breaks=c(10, 100, 500, 1000, 10000, 50000)) +
  ylab("Millions") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

r.fig



shp.mp <-readOGR(dsn="~/Projects/GOA ichthyoplankton time series/NewLateLarvalPolygon",layer="Survey_poly")
shp.mp.LL<-spTransform(shp.mp,CRS("+proj=longlat"))

map('worldHires',xlim=c(-170,-140),ylim=c(52,62),fill=T,col="gray",border=F)
map(shp.mp.LL,add=T,col="blue",lwd=2)
map.axes()

### Late summer age-0 trawl POLYGON: doesn't look quite as tidy, but this is what I use to extract the trawl stations for the "core" area from our late summer survey. 
### I've also included a pdf showing the last 4 years of pollock catches. But in the 2000s, stations were restricted mostly to the "core" area, hence our subsetting for the timeseries.
polygon(c(-159.7,-158.3,-155,-156.3),c(55.9,54.7,56.1,57.2))