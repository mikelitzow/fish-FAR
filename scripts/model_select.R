## Model selection tables

library(rstan)

## 1. ZINB cod seine models (fig. 2a) ----------------------
## y = cod CPUE
cod0_zinb_k3   <- readRDS("./output/cod0_zinb_k3.rds")   ## bay_fac
cod0s_zinb_k3  <- readRDS("./output/cod0s_zinb_k3.rds")  ## bay_fac/site_fac
cod1sg_zinb_k3 <- readRDS("./output/cod1sg_zinb_k3.rds") ## bay_fac/site_fac + temp
cod2sg_zinb_k3 <- readRDS("./output/cod2sg_zinb_k3.rds") ## bay_fac/site_fac + temp + ssb

loo(cod0_zinb_k3, cod0s_zinb_k3,
    cod1sg_zinb_k3, cod2sg_zinb_k3)

m1_looic    <- c(cod0_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod0s_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod1sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod2sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"])
m1_name     <- rep("M1", length(m1_looic))
m1_response <- rep("Cod CPUE", length(m1_looic))
m1_covars   <- c("bay", "bay/site", "bay/site + temp", "bay/site + temp + ssb")
m1_delta    <- m1_looic - min(m1_looic)
m1_tab      <- data.frame(name = m1_name,
                     response = m1_response,
                     covars = m1_covars,
                     LOOIC = m1_looic,
                     LOOIC_delta = m1_delta)


## 2. ZINB cod seine models (fig 2c.) ----------------------
## y = cod CPUE
recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")  ## year_fac + bay_fac
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")  ## year_fac + bay_fac/site_fac

loo(recr_1_zinb, recr_2_zinb)

m2_looic    <- c(recr_1_zinb$criteria$loo$estimates["looic", "Estimate"],
              recr_2_zinb$criteria$loo$estimates["looic", "Estimate"])
m2_name     <- rep("M2", length(m2_looic))
m2_response <- rep("Cod CPUE", length(m2_looic))
m2_covars   <- c("year + bay", "year + bay/site")
m2_delta    <- m2_looic - min(m2_looic)
m2_tab      <- data.frame(name = m2_name,
                     response = m2_response,
                     covars = m2_covars,
                     LOOIC = m2_looic,
                     LOOIC_delta = m2_delta)


## 3. Cod recruitment models (fig. 2c) ---------------------
## (predicting stock assessment model recruitment from seine data)
## y = assessment recruitment
codR1_brm  <- readRDS("./output/codR1_brm.rds") ## seine + ssb
codR2_brm  <- readRDS("./output/codR2_brm.rds") ## seine

loo(codR1_brm, codR2_brm)

m3_looic    <- c(codR1_brm$criteria$loo$estimates["looic", "Estimate"],
              codR2_brm$criteria$loo$estimates["looic", "Estimate"])
m3_name     <- rep("M3", length(m3_looic))
m3_response <- rep("Recruitment", length(m3_looic))
m3_covars   <- c("seine + ssb", "seine")
m3_delta    <- m3_looic - min(m3_looic)
m3_tab      <- data.frame(name = m3_name,
                     response = m3_response,
                     covars = m3_covars,
                     LOOIC = m3_looic,
                     LOOIC_delta = m3_delta)

## 4. Pollock DFA (fig. 3a) --------------------------------
## y = pollock DFA trend
dfa1_far_brm  <- readRDS("./output/dfa1_far_brm.rds")  ## ssb + far
dfa2_far_brm  <- readRDS("./output/dfa2_far_brm.rds")  ## far

loo(dfa1_far_brm, dfa2_far_brm)

m4_looic    <- c(dfa1_far_brm$criteria$loo$estimates["looic", "Estimate"],
              dfa2_far_brm$criteria$loo$estimates["looic", "Estimate"])
m4_name     <- rep("M4", length(m4_looic))
m4_response <- rep("DFA trend", length(m4_looic))
m4_covars   <- c("far + ssb", "far")
m4_delta    <- m4_looic - min(m4_looic)
m4_tab      <- data.frame(name = m4_name,
                     response = m4_response,
                     covars = m4_covars,
                     LOOIC = m4_looic,
                     LOOIC_delta = m4_delta)

## dfa_temp1_brm  <- readRDS("./output/dfa_temp1_brm.rds") ## ssb + mean.anom
## dfa_temp2_brm  <- readRDS("./output/dfa_temp2_brm.rds") ## mean.anom
## loo(dfa_temp1_brm, dfa_temp2_brm)


## 5. Pollock recruit projections (Fig. 2c_3b_4b) ----------
## y = pollock recruitment
poll.R1 <- readRDS("./output/poll_R1_FAR_obs.rds")   ## far
poll.R1s <- readRDS("./output/poll_R1s_FAR_obs.rds") ## far + ssb

loo(poll.R1, poll.R1s)

m5_looic    <- c(poll.R1$criteria$loo$estimates["looic", "Estimate"],
              poll.R1s$criteria$loo$estimates["looic", "Estimate"])
m5_name     <- rep("M5", length(m5_looic))
m5_response <- rep("Recruitment", length(m5_looic))
m5_covars   <- c("far", "far + ssb")
m5_delta    <- m5_looic - min(m5_looic)
m5_tab      <- data.frame(name = m5_name,
                     response = m5_response,
                     covars = m5_covars,
                     LOOIC = m5_looic,
                     LOOIC_delta = m5_delta)


## Final Table ---------------------------------------------
model_table <- rbind(m1_tab, m2_tab, m3_tab, m4_tab, m5_tab)
write.csv(model_table, "./output/model_table.csv", row.names = FALSE)
