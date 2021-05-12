## Model selection tables

library(rstan)
library(brms)

## 1. ZINB cod seine models (fig. 2a) ----------------------
## y = cod CPUE
cod0_zinb_k3   <- readRDS("./output/cod0_zinb_k3.rds")   ## bay_fac
cod0s_zinb_k3  <- readRDS("./output/cod0s_zinb_k3.rds")  ## bay_fac/site_fac
cod1sg_zinb_k3 <- readRDS("./output/cod1sg_zinb_k3.rds") ## bay_fac/site_fac + temp
cod2sg_zinb_k3 <- readRDS("./output/cod2sg_zinb_k3.rds") ## bay_fac/site_fac + temp + SSB

loo(cod0_zinb_k3, cod0s_zinb_k3,
    cod1sg_zinb_k3, cod2sg_zinb_k3)

m1_looic    <- c(cod0_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod0s_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod1sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"],
              cod2sg_zinb_k3$criteria$loo$estimates["looic", "Estimate"])
m1_name     <- rep("M1", length(m1_looic))
m1_response <- rep("Cod seine abundance", length(m1_looic))
m1_covars   <- c("DOY + bay",
                 "DOY + bay/site",
                 "DOY + bay/site + temp",
                 "DOY + bay/site + temp + SSB")
m1_delta    <- m1_looic - min(m1_looic)
m1_tab      <- data.frame(name = m1_name,
                          response = m1_response,
                          note = "Fig. 2a",
                          covars = m1_covars,
                          LOOIC = m1_looic,
                          LOOIC_delta = m1_delta)
m1_tab <- m1_tab[order(m1_tab$LOOIC_delta), ]
m1_tab$response <- c(m1_tab$response[1], rep("", nrow(m1_tab) - 1))
m1_tab$note <- c(m1_tab$note[1], rep("", nrow(m1_tab) - 1))

## Best model
m1_mat <- as.matrix(cod2sg_zinb_k3, pars = c("^Intercept", "sds_", "sd"))
m1_sum <- t(apply(m1_mat, 2, function(x) quantile(x, probs = c(0.0275, 0.5, 0.975))))
m1_best <- cbind(name = "M1", parameter = row.names(m1_sum), stringsAsFactors = F, as.data.frame(m1_sum))
row.names(m1_best) <- NULL
m1_pars <- list(c("Intercept", "Intercept"),
                c("Intercept_zi", "ZI Intercept"),
                c("sds_sjulian_1", "SD Smooth: DOY"),
                c("sds_stemp.anom_1", "SD Smooth: Temp"),
                c("sds_sssb_1", "SD Smooth: SSB"),
                c("sds_zi_sjulian_1", "ZI SD Smooth: DOY"),
                c("sds_zi_stemp.anom_1", "ZI SD Smooth: Temp"),
                c("sds_zi_sssb_1", "ZI SD Smooth: SSB"),
                c("sd_bay_fac__Intercept", "SD: bay effect"),
                c("sd_bay_fac:site_fac__Intercept", "SD: site effect"),
                c("sd_bay_fac__zi_Intercept", "ZI SD: bay effect"),
                c("sd_bay_fac:site_fac__zi_Intercept", "ZI SD: site effect"))
for(i in seq_along(m1_pars)) {
    m1_best$parameter[m1_best$parameter == m1_pars[[i]][1]] <- m1_pars[[i]][2]
}


## 2. ZINB cod seine models (fig 2c.) ----------------------
## y = cod CPUE
recr_1_zinb <- readRDS("./output/recr_1_zinb.rds")  ## year_fac + bay_fac
recr_2_zinb <- readRDS("./output/recr_2_zinb.rds")  ## year_fac + bay_fac/site_fac

loo(recr_1_zinb, recr_2_zinb)

m2_looic    <- c(recr_1_zinb$criteria$loo$estimates["looic", "Estimate"],
              recr_2_zinb$criteria$loo$estimates["looic", "Estimate"])
m2_name     <- rep("M2", length(m2_looic))
m2_response <- rep("Cod annual seine abundance", length(m2_looic))
m2_covars   <- c("year + DOY + bay", "year + DOY + bay/site")
m2_delta    <- m2_looic - min(m2_looic)
m2_tab      <- data.frame(name = m2_name,
                          response = m2_response,
                          note = "For comparing to stock assessment model",
                          covars = m2_covars,
                          LOOIC = m2_looic,
                          LOOIC_delta = m2_delta)
m2_tab <- m2_tab[order(m2_tab$LOOIC_delta), ]
m2_tab$response <- c(m2_tab$response[1], rep("", nrow(m2_tab) - 1))
m2_tab$note <- c(m2_tab$note[1], rep("", nrow(m2_tab) - 1))

## Best model
m2_mat <- as.matrix(recr_2_zinb, pars = c("^Intercept", "sds_", "sd"))
m2_sum <- t(apply(m2_mat, 2, function(x) quantile(x, probs = c(0.0275, 0.5, 0.975))))
m2_best <- cbind(name = "M2", parameter = row.names(m2_sum), stringsAsFactors = F, as.data.frame(m2_sum))
row.names(m2_best) <- NULL
m2_pars <- list(c("Intercept", "Intercept"),
                c("Intercept_zi", "ZI Intercept"),
                c("sds_sjulian_1", "SD Smooth: DOY"),
                c("sds_zi_sjulian_1", "ZI SD Smooth: DOY"),
                c("sd_bay_fac__Intercept", "SD: bay effect"),
                c("sd_bay_fac:site_fac__Intercept", "SD: site effect"),
                c("sd_bay_fac__zi_Intercept", "ZI SD: bay effect"),
                c("sd_bay_fac:site_fac__zi_Intercept", "ZI SD: site effect"))
for(i in seq_along(m2_pars)) {
    m2_best$parameter[m2_best$parameter == m2_pars[[i]][1]] <- m2_pars[[i]][2]
}


## 3. Cod recruitment models (fig. 2c) ---------------------
## (predicting stock assessment model recruitment from seine data)
## y = assessment recruitment
codR1_brm  <- readRDS("./output/codR1_brm.rds") ## seine + SSB
codR2_brm  <- readRDS("./output/codR2_brm.rds") ## seine

loo(codR1_brm, codR2_brm)

m3_looic    <- c(codR1_brm$criteria$loo$estimates["looic", "Estimate"],
              codR2_brm$criteria$loo$estimates["looic", "Estimate"])
m3_name     <- rep("M3", length(m3_looic))
m3_response <- rep("Cod stock assessment model recruitment", length(m3_looic))
m3_covars   <- c("seine + SSB", "seine")
m3_delta    <- m3_looic - min(m3_looic)
m3_tab      <- data.frame(name = m3_name,
                          response = m3_response,
                          note = "To estimate 2017-2020 values in Fig. 2c",
                          covars = m3_covars,
                          LOOIC = m3_looic,
                          LOOIC_delta = m3_delta)
m3_tab <- m3_tab[order(m3_tab$LOOIC_delta), ]
m3_tab$response <- c(m3_tab$response[1], rep("", nrow(m3_tab) - 1))
m3_tab$note <- c(m3_tab$note[1], rep("", nrow(m3_tab) - 1))

## Best model
m3_mat <- as.matrix(codR2_brm, pars = c("^Intercept", "seine"))
m3_sum <- t(apply(m3_mat, 2, function(x) quantile(x, probs = c(0.0275, 0.5, 0.975))))
m3_best <- cbind(name = "M3", parameter = row.names(m3_sum), stringsAsFactors = F, as.data.frame(m3_sum))
row.names(m3_best) <- NULL
m3_pars <- list(c("Intercept", "Intercept"),
                c("b_seine", "Seine"))
for(i in seq_along(m3_pars)) {
    m3_best$parameter[m3_best$parameter == m3_pars[[i]][1]] <- m3_pars[[i]][2]
}



## 4. Pollock DFA (fig. 3a) --------------------------------
## y = pollock DFA trend
dfa1_far_brm  <- readRDS("./output/dfa1_far_brm.rds")  ## SSB + far
dfa2_far_brm  <- readRDS("./output/dfa2_far_brm.rds")  ## far
dfa_temp1_brm  <- readRDS("./output/dfa_temp1_brm.rds") ## SSB + mean.anom
dfa_temp2_brm  <- readRDS("./output/dfa_temp2_brm.rds") ## mean.anom
dfa_temp3_brm  <- readRDS("./output/dfa_temp3_brm.rds") ## SSB + larval anom
dfa_temp4_brm  <- readRDS("./output/dfa_temp4_brm.rds") ## larval.anom

loo(dfa1_far_brm, dfa2_far_brm, dfa_temp1_brm, dfa_temp2_brm, dfa_temp3_brm, dfa_temp4_brm)


m4_looic    <- c(dfa1_far_brm$criteria$loo$estimates["looic", "Estimate"],
                 dfa2_far_brm$criteria$loo$estimates["looic", "Estimate"],
                 dfa_temp1_brm$criteria$loo$estimates["looic", "Estimate"],
                 dfa_temp2_brm$criteria$loo$estimates["looic", "Estimate"],
                 dfa_temp3_brm$criteria$loo$estimates["looic", "Estimate"],
                 dfa_temp4_brm$criteria$loo$estimates["looic", "Estimate"])

m4_name     <- rep("M4", length(m4_looic))
m4_response <- rep("Pollock DFA trend", length(m4_looic))
m4_covars   <- c("FAR + SSB", "FAR",
                 "SSB + egg/larval temp", "egg/larval temp",
                 "SSB + larval temp", "larval temp")
m4_delta    <- m4_looic - min(m4_looic)
m4_tab      <- data.frame(name = m4_name,
                          response = m4_response,
                          note = "Fig. 3a",
                          covars = m4_covars,
                          LOOIC = m4_looic,
                          LOOIC_delta = m4_delta)
m4_tab <- m4_tab[order(m4_tab$LOOIC_delta), ]
m4_tab$response <- c(m4_tab$response[1], rep("", nrow(m4_tab) - 1))
m4_tab$note <- c(m4_tab$note[1], rep("", nrow(m4_tab) - 1))

## Best model
m4_mat <- as.matrix(dfa2_far_brm, pars = c("^Intercept", "sds_", "sd"))
m4_sum <- t(apply(m4_mat, 2, function(x) quantile(x, probs = c(0.0275, 0.5, 0.975))))
m4_best <- cbind(name = "M4", parameter = row.names(m4_sum), stringsAsFactors = F, as.data.frame(m4_sum))
row.names(m4_best) <- NULL
m4_pars <- list(c("Intercept", "Intercept"),
                c("sds_sfar_1", "SD Smooth: FAR"))
for(i in seq_along(m4_pars)) {
    m4_best$parameter[m4_best$parameter == m4_pars[[i]][1]] <- m4_pars[[i]][2]
}



## 5. Pollock recruit projections (Fig. 2c_3b_4b) ----------
## y = pollock recruitment
poll.R1 <- readRDS("./output/poll_R1_FAR_obs.rds")   ## far
poll.R1s <- readRDS("./output/poll_R1s_FAR_obs.rds") ## far + SSB

loo(poll.R1, poll.R1s)

m5_looic    <- c(poll.R1$criteria$loo$estimates["looic", "Estimate"],
              poll.R1s$criteria$loo$estimates["looic", "Estimate"])
m5_name     <- rep("M5", length(m5_looic))
m5_response <- rep("Pollock model recruitment estimate", length(m5_looic))
m5_covars   <- c("FAR", "FAR + SSB")
m5_delta    <- m5_looic - min(m5_looic)
m5_tab      <- data.frame(name = m5_name,
                          response = m5_response,
                          note = "Fig. 3b",
                          covars = m5_covars,
                          LOOIC = m5_looic,
                          LOOIC_delta = m5_delta)
m5_tab <- m5_tab[order(m5_tab$LOOIC_delta), ]
m5_tab$response <- c(m5_tab$response[1], rep("", nrow(m5_tab) - 1))
m5_tab$note <- c(m5_tab$note[1], rep("", nrow(m5_tab) - 1))

## Best model
m5_mat <- as.matrix(poll.R1, pars = c("^Intercept", "sds_", "sd"))
m5_sum <- t(apply(m5_mat, 2, function(x) quantile(x, probs = c(0.0275, 0.5, 0.975))))
m5_best <- cbind(name = "M5", parameter = row.names(m5_sum), stringsAsFactors = F, as.data.frame(m5_sum))
row.names(m5_best) <- NULL
m5_pars <- list(c("Intercept", "Intercept"),
                c("sds_sFAR_1", "SD Smooth: FAR"))
for(i in seq_along(m5_pars)) {
    m5_best$parameter[m5_best$parameter == m5_pars[[i]][1]] <- m5_pars[[i]][2]
}


## Final Table ---------------------------------------------
model_select <- rbind(m1_tab, m2_tab, m3_tab, m4_tab, m5_tab)
names(model_select) <- c("Model set", "Response variable", "Notes",
                         "Covariates", "LOOIC", "delta_LOOIC")
write.csv(model_select, "./output/model_select.csv", row.names = FALSE)

model_best <- rbind(m1_best, m2_best, m3_best, m4_best, m5_best)
names(model_best) <- c("Model set", "Parameter", "2.5%", "50%", "97.5%")
write.csv(model_best, "./output/model_best.csv", row.names = FALSE)
