library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(brms)
library(rgdal)

theme_set(theme_bw())

# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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

far.fig <- ggplot(filter(pred.obs, year >= 1970)) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90") +
  geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  theme(axis.title.x = element_blank()) +
  ylab("FAR") +
  scale_x_continuous(breaks=seq(1960, 2020, 10)) 

print(far.fig)


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

## study site -----------------------------------------------

# load data sets
bays <- read.csv("./data/bay_lat_long.csv", row.names = 1)

# remove Cooks (repeat) and Kujulik (not resampled)
drop <- bays$Bay %in% c("Cooks", "Kujulik")
bays <- bays[!drop,]
bays$years <- ifelse(bays$Bay %in% c("Anton Larson Bay", "Cook Bay"), "2006-2020", "2018-2020")


ak <- ne_countries(scale = "large", returnclass = "sf", continent="north america")

# use this version unless the high-res version is entered!
# ak <- ne_countries(scale = "medium", returnclass = "sf", continent="north america")
world <- ne_countries(scale='medium', returnclass = "sf")

# add FOCI
shp.mp <-readOGR(dsn="./data/FOCI_survey_polygon",layer="Survey_poly")
shp.mp.LL<-spTransform(shp.mp,CRS("+proj=longlat"))


## need to change Spatial Polygon to dataframe
# add to data a new column termed "id" composed of the rownames of data
shp.mp.LL@data$id <- rownames(shp.mp.LL@data)

# create a data.frame from our spatial object
poly.points <- fortify(shp.mp.LL, region = "id")

# merge the "fortified" data with the data from our spatial object
ichthyo <- merge(poly.points, shp.mp.LL@data, by = "id")

# combine with age-0 trawl polygon
trawl <- data.frame(long = c(-159.7,-158.3,-155,-156.3,-159.7),
                    lat = c(55.9,54.7,56.1,57.2,55.9),
                    type = "Juvenile trawl")

ichthyo <- ichthyo %>%
  select(long, lat) %>%
  mutate(type = "Larval survey")

polys <- rbind(ichthyo, trawl)
polys$type <- reorder(polys$type, desc(polys$type))

bays$type <- "Beach seine"


box <- data.frame(long = c(-163, -163, -151, -151, -163), lat = c(54.5, 58.5, 58.5, 54.5, 54.5))

inset <- ggplot(data = world) +
  geom_sf(fill="dark grey", color=NA) +
  coord_sf(xlim = c(-179, -70), ylim = c(0, 70)) +
  geom_path(data=box, aes(long, lat), color=cb[8], size=1) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(color="black", fill="transparent"),
        panel.spacing = unit(1, 'mm'))

inset  



map.plot <- ggplot(ak) +  
  geom_path(data=polys, aes(long, lat, color=type), lwd=1.5) +
  geom_sf(fill="darkgoldenrod3", color=NA) + 
  coord_sf(xlim = c(-163, -151), ylim = c(54.5, 58.5), expand = FALSE) +
  geom_point(data = bays, aes(-lon, lat, fill=type), size=3, shape=21) +
  theme(axis.title = element_blank(),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.margin = margin(-2,0,0,0,unit='mm'),
        legend.background = element_rect(fill = 'transparent', linetype=0),
        legend.spacing.y = unit(1, 'mm')) +
  scale_fill_manual(values=cb[4]) +
  scale_color_manual(values=cb[c(6,7)]) +
  scale_x_continuous(breaks = c(-160, -156, -152)) +
  scale_y_continuous(breaks = c(55, 56, 57, 58))

map.plot

full.map <- map.plot +
  annotation_custom(
    grob = ggplotGrob(inset),
    xmin = -163,
    xmax = -159,
    ymin = 56,
    ymax = 59
  ) 

full.map

png("./figs/Fig1.png", width = 4, height = 7, units = 'in', res = 300)
ggpubr::ggarrange(far.fig, ssb.fig, full.map, 
                  ncol=1, 
                  labels = c("a", "b", "c"))
dev.off()
