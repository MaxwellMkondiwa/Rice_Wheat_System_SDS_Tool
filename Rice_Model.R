# Rice Model

library(sp)
library(mgcv)
library(bamlss)

# Multivariate geoadditive model
# remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")
#
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

f_rice_yield_MRF <- list(
  b_grain_yield_ton_per_ha_rice ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  sigma ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re")
)

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
b_rice_yield_MRF <- bamlss(f_rice_yield_MRF, data = Irrig_Rev_rice_wheat, family = "gaussian")

library(distreg.vis)


if (interactive()) {
   distreg.vis:: vis()
 }