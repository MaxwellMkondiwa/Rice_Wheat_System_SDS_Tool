

#rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(sp, mgcv, bamlss,BayesX,R2BayesX,sf,spdep,rio,distreg.vis)

setwd("C:/Users/MMKONDIWA/OneDrive - CIMMYT/Documents/GitHub/Rice_Wheat_System_SDS_Tool")
# Multivariate geoadditive model
# remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")
#
library(sp)
library (mgcv)
library(mvnchol)
library(BayesX)
library(R2BayesX)
library(sf)
library(spdep)
library(rio)

Irrig_Rev_rice_wheat <- import("data/Irrig_Rev_rice_wheat_Updated.csv")

Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$Longitude)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$Latitude)))

library(janitor)
library(tidyr)


Irrig_Rev_rice_wheat$harvest_day_rice <- Irrig_Rev_rice_wheat$l_crop_duration_days_rice + Irrig_Rev_rice_wheat$sowdate_fmt_rice_day

shpname <- file.path(getwd(), "shp", "India_aoi_sf_sp")

India_aoi_sp_bnd <- BayesX::shp2bnd(shpname = shpname, regionnames = "District", check.is.in = F)
# Multivariate geoadditive model
# remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")
# install_github("https://github.com/meteosimon/mvnchol")
library(mvnchol)

# Remove NAs in the monsoon variables
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$onset_2017)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$monsoon_onset_dev)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$median_onset_82_15)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$sd_onset_82_15)))


K <- neighbormatrix(India_aoi_sp_bnd)
head(K)

## Also need to transform to factor for
## setting up the MRF smooth.
Irrig_Rev_rice_wheat$District <- as.factor(Irrig_Rev_rice_wheat$a_q103_district)

## Now note that not all regions are observed,
## therefore we need to remove those regions
## from the penalty matrix
rn <- rownames(K)
lv <- levels(Irrig_Rev_rice_wheat$District)
i <- rn %in% lv
K <- K[i, i]

set.seed(321)
library(bamlss)


# Rice first stage
f_sow_rice_1st_stage <- list(
  sowdate_fmt_rice_day ~ 1 + rice_duration_class_long + s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re")
)

f_sow_rice_1st_stage_MRF <- bamlss(f_sow_rice_1st_stage, data = Irrig_Rev_rice_wheat, family = "gaussian")

fitted_f_sow_rice_1st_stage_MRF <- f_sow_rice_1st_stage_MRF$fitted

Irrig_Rev_rice_wheat$Res_rice_sow <- Irrig_Rev_rice_wheat$sowdate_fmt_rice_day - fitted_f_sow_rice_1st_stage_MRF$mu

# Wheat first stage
f_sow_wheat_1st_stage <- list(
  sowdate_fmt_wheat_day ~ 1 + variety_type_NMWV + s(harvest_day_rice) + s(gw_2018) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re")
)

f_sow_wheat_1st_stage_MRF <- bamlss(f_sow_wheat_1st_stage, data = Irrig_Rev_rice_wheat, family = "gaussian")

fitted_f_sow_wheat_1st_stage_MRF <- f_sow_wheat_1st_stage_MRF$fitted

Irrig_Rev_rice_wheat$Res_wheat_sow <- Irrig_Rev_rice_wheat$sowdate_fmt_wheat_day - fitted_f_sow_wheat_1st_stage_MRF$mu


# Multivariate with sowing dates as endogenous
f_rice_wheat_yield_MRF_corr_endo <- list(
  sowdate_fmt_rice_day ~ 1 + rice_duration_class_long + s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  b_grain_yield_ton_per_ha_rice ~ 1 + rice_duration_class_long + s(Res_rice_sow) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  sowdate_fmt_wheat_day ~ 1 + variety_type_NMWV + s(harvest_day_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  l_ton_per_hectare ~ 1 + variety_type_NMWV + s(Res_wheat_sow) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lamdiag1 ~ 1,
  lamdiag2 ~ 1,
  lamdiag3 ~ 1,
  lamdiag4 ~ 1,
  lambda12 ~ 1,
  lambda13 ~ 1,
  lambda14 ~ 1,
  lambda23 ~ 1,
  lambda24 ~ s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lambda34 ~ 1
)

multivariate_geo_sow_MRF_corr_endo <- bamlss(f_rice_wheat_yield_MRF_corr_endo, type = "modified", family = mvnchol_bamlss(k = 4), data = Irrig_Rev_rice_wheat)

library(distreg.vis)
library(bamlss)

library(mvnchol)
if (interactive()) {
  distreg.vis:: vis()
}