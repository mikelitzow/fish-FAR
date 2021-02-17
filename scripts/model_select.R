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


## 2. ZINB cod seine models (fig 2c.) ----------------------
## y = cod CPUE
recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")  ## year_fac + bay_fac
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")  ## year_fac + bay_fac/site_fac

loo(recr_1_zinb, recr_2_zinb)


## 3. Cod recruitment models (fig. 2c) ---------------------
## (predicting stock assessment model recruitment from seine data)
## y = assessment recruitment
codR1_brm  <- readRDS("./output/codR1_brm.rds") ## seine + ssb
codR2_brm  <- readRDS("./output/codR2_brm.rds") ## seine

loo(codR1_brm, codR2_brm)


## 4. Pollock DFA (fig. 3a) --------------------------------
## y = pollock DFA trend
dfa1_far_brm  <- readRDS("./output/dfa1_far_brm.rds")  ## ssb + far
dfa2_far_brm  <- readRDS("./output/dfa2_far_brm.rds")  ## far

loo(dfa1_far_brm, dfa2_far_brm)

dfa_temp1_brm  <- readRDS("./output/dfa_temp1_brm.rds") ## ssb + mean.anom
dfa_temp2_brm  <- readRDS("./output/dfa_temp2_brm.rds") ## mean.anom

loo(dfa_temp1_brm, dfa_temp2_brm)


## 5. Pollock recruit projections (Fig. 2c_3b_4b) ----------
## y = pollock recruitment
poll.R1 <- readRDS("./output/poll_R1_FAR_obs.rds")   ## far
poll.R1s <- readRDS("./output/poll_R1s_FAR_obs.rds") ## far + ssb

loo(poll.R1, poll.R1s)
