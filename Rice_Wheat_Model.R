

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



# Yield analysis

## Non structural ------------------------------------
f_rice_wheat_yield_MRF_corr_endo_yld_non_struct<- list(
  b_grain_yield_ton_per_ha_rice ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  l_ton_per_hectare ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lamdiag1 ~ 1,
  lamdiag2 ~ 1,
  lambda12 ~ 1+s(sowdate_fmt_rice_day)+s(sowdate_fmt_wheat_day)+g_q5305_irrig_times_rice+rice_duration_class_long+variety_type_NMWV+s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re")
)


multivariate_geo_sow_MRF_corr_endo_yld_non_struct <- bamlss(f_rice_wheat_yield_MRF_corr_endo_yld_non_struct, type = "modified", family = mvnchol_bamlss(k = 2), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo_yld_non_struct)

plot(multivariate_geo_sow_MRF_corr_endo_yld_non_struct)



## Structural --------------------------------------

### First stage ------------------------------------------
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


# Multivariate with sowing dates as endogenous ---------------------
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
  lambda24 ~ 1+s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lambda34 ~ 1
)

multivariate_geo_sow_MRF_corr_endo <- bamlss(f_rice_wheat_yield_MRF_corr_endo, type = "modified", family = mvnchol_bamlss(k = 4), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo)

plot(multivariate_geo_sow_MRF_corr_endo)


## More results -------------------------------
nd <- data.frame("District" = unique(Irrig_Rev_rice_wheat$District))

# Focusing on the cross-equation correlations

## Predict for the structured spatial effects.
p_str_multivariate_geo_sow_MRF_corr_endo_rice_y <- predict(multivariate_geo_sow_MRF_corr_endo, newdata = nd, term = "s(District,id='mrf1')", intercept = FALSE)
p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt <- as.data.frame(p_str_multivariate_geo_sow_MRF_corr_endo_rice_y)
write.csv(p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt, "tables/p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt2.csv")

## And the unstructured spatial effect.
p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y <- predict(multivariate_geo_sow_MRF_corr_endo, newdata = nd, term = "s(District,id='re2')", intercept = FALSE)

# Rice sowing spatial equation

plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y$mu1, id = nd$District, main = expression(mu(rice_sowing)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y$mu1, id = nd$District, main = expression(mu(rice_sowing)), title = "Unstructured spatial effect")

# Rice yield spatial equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y$mu2, id = nd$District, main = expression(mu(rice_yield)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y$mu2, id = nd$District, main = expression(mu(rice_yield)), title = "Unstructured spatial effect")

# Wheat sowing equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y$mu3, id = nd$District, main = expression(mu(wheat_sowing)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y$mu3, id = nd$District, main = expression(mu(wheat_sowing)), title = "Unstructured spatial effect")

# Wheat yield spatial equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y$mu4, id = nd$District, main = expression(mu(wheat_yield)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y$mu4, id = nd$District, main = expression(mu(wheat_yield)), title = "Unstructured spatial effect")

# Focusing on the cross-equation correlations

## MRF smooth effect.
plotmap(India_aoi_sp_bnd,
        x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y$lambda24, id = nd$District,
        main = expression(lambda(rice, wheat)), title = "Structured spatial effect"
)

## Random effects.
plotmap(India_aoi_sp_bnd,
        x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y$lambda24, id = nd$District, main = expression(lambda(rice, wheat)), title = "Unstructured spatial effect"
)

# Rice sowing equation : Non-linear relationships
# s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15)
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu1", term = "s(gw_2017)")
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu1", term = "s(onset_2017)")
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu1", term = "s(monsoon_onset_dev)")
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu1", term = "s(median_onset_82_15)")
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu1", term = "s(sd_onset_82_15)")


# Rice yield equation
# s(Res_rice_sow) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice)

plot(multivariate_geo_sow_MRF_corr_endo, model = "mu2", term = "s(Res_rice_sow)")

plot(multivariate_geo_sow_MRF_corr_endo, model = "mu2", term = "s(g_q5305_irrig_times_rice)")

plot(multivariate_geo_sow_MRF_corr_endo, model = "mu2", term = "s(nperha_rice)")

plot(multivariate_geo_sow_MRF_corr_endo, model = "mu2", term = "s(p2o5perha_rice)")

# Wheat sowing equation
# s(harvest_day_rice) + s(gw_2018)
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu3", term = "s(harvest_day_rice)")

# (multivariate_geo_sow_MRF_corr_endo, model = "mu3", term = "s(gw_2018)")

# Wheat yield equation
plot(multivariate_geo_sow_MRF_corr_endo, model = "mu4", term = "s(Res_wheat_sow)")


# Fitted values
multivariate_geo_sow_MRF_corr_endo_fitted_values <- multivariate_geo_sow_MRF_corr_endo$fitted

multivariate_geo_sow_MRF_corr_endo_fitted_values <- as.data.frame(multivariate_geo_sow_MRF_corr_endo_fitted_values)

# Merge the fitted results to the data and export

Irrig_Rev_rice_wheat_Mult_Res <- cbind(Irrig_Rev_rice_wheat, multivariate_geo_sow_MRF_corr_endo_fitted_values)

write.csv(Irrig_Rev_rice_wheat_Mult_Res, "tables/Irrig_Rev_rice_wheat_Mult_Res.csv")

summary(lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res))

summary(lm(b_grain_yield_ton_per_ha_rice ~ l_ton_per_hectare, Irrig_Rev_rice_wheat_Mult_Res))



plot((lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res)))

Irrig_Rev_rice_wheat_Mult_Res_sp <- Irrig_Rev_rice_wheat_Mult_Res
coordinates(Irrig_Rev_rice_wheat_Mult_Res_sp) <- c("o_largest_plot_gps_longitude", "o_largest_plot_gps_latitude")

# Map the correlations
library(modelsummary)
Irrig_Rev_rice_wheat_Mult_Res_dist <- datasummary(Heading("District") * District ~ Heading("N obs") * N + Heading("%") * Percent() + lambda24 * (Mean + SD), data = Irrig_Rev_rice_wheat_Mult_Res, output = "data.frame")

library(dplyr)
Irrig_Rev_rice_wheat_Mult_Res_dist <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist, "Mean_Rice_Wheat_Rho" = "Mean")
Irrig_Rev_rice_wheat_Mult_Res_dist <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist, "SD_Rice_Wheat_Rho" = "SD")

Irrig_Rev_rice_wheat_Mult_Res_dist <- subset(Irrig_Rev_rice_wheat_Mult_Res_dist, Irrig_Rev_rice_wheat_Mult_Res_dist$District != "Purnia")


# Bihar and EUP map
# India district Map

library(geodata)

India <- gadm(country = "IND", level = 2, path = "shp")

plot(India)

India_aoi <- subset(India, India$NAME_1 == "Bihar" | India$NAME_2 %in% c("Ballia", "Chandauli", "Deoria", "Ghazipur", "Kushinagar", "Maharajganj", "Mau", "Siddharth Nagar", "Gorakhpur"))

plot(India_aoi)

plot(India_aoi, add = TRUE)

library(sf)

India_aoi_sf <- st_as_sf(India_aoi)
library(mapview)

mapview(India_aoi_sf)

# Dissolve the district polygons to form new polygon of Bihar and EUP
library(sf)
India_aoi_sf_dis <- st_union(India_aoi_sf)
mapview(India_aoi_sf_dis)

###
library(modelsummary)

India_aoi_sf$District <- India_aoi_sf$NAME_2

India_aoi_sf$District[India_aoi_sf$District == "Purba Champaran"] <- "EastChamparan"
India_aoi_sf$District[India_aoi_sf$District == "Pashchim Champaran"] <- "WestChamparan"
India_aoi_sf$District[India_aoi_sf$District == "Bhojpur"] <- "Arah"
India_aoi_sf$District[India_aoi_sf$District == "Ballia"] <- "Balia"
India_aoi_sf$District[India_aoi_sf$District == "Ghazipur"] <- "Gazipur"
India_aoi_sf$District[India_aoi_sf$District == "Siddharth Nagar"] <- "Siddharthnagar"

rice_wheat_yield_rho_dist_sf <- merge(India_aoi_sf, Irrig_Rev_rice_wheat_Mult_Res_dist, by = "District")

rice_wheat_yield_rho_dist_sf$Mean_Rice_Wheat_Rho <- as.numeric(rice_wheat_yield_rho_dist_sf$Mean_Rice_Wheat_Rho)
library(mapview)
mapview(rice_wheat_yield_rho_dist_sf, zcol = "Mean_Rice_Wheat_Rho", layer.name = "Rice wheat equation correlation")
library(sf)
rice_wheat_yield_rho_dist_sf_sp <- as_Spatial(rice_wheat_yield_rho_dist_sf)
library(tmap)
tmap_mode("plot")
rice_wheat_yield_rho_dist_sf_sp_map <- tm_shape(rice_wheat_yield_rho_dist_sf_sp) +
  tm_polygons(col = "Mean_Rice_Wheat_Rho", title = "Rice wheat equation correlation",style="quantile",palette = "YlGn") +
  tm_layout(legend.outside = TRUE)

rice_wheat_yield_rho_dist_sf_sp_map

tmap_save(rice_wheat_yield_rho_dist_sf_sp_map, "figures/rice_wheat_yield_rho_dist_sf_sp_map .png")





# Revenue analysis -----------------------------------------------
### Non structural ---------------------------------
f_rice_wheat_yield_MRF_corr_endo_rev_non_struct<- list(
  revenue_rice ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  revenue_wheat ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lamdiag1 ~ 1,
  lamdiag2 ~ 1,
  lambda12 ~ 1+s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re")
)


multivariate_geo_sow_MRF_corr_endo_rev_non_struct <- bamlss(f_rice_wheat_yield_MRF_corr_endo_rev_non_struct, type = "modified", family = mvnchol_bamlss(k = 2), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo_rev_non_struct)

plot(multivariate_geo_sow_MRF_corr_endo_rev_non_struct)






### Structural -----------------------------------
f_rice_wheat_yield_MRF_corr_endo_rev <- list(
  sowdate_fmt_rice_day ~ 1 + rice_duration_class_long + s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  revenue_rice ~ 1 + rice_duration_class_long + s(Res_rice_sow) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  sowdate_fmt_wheat_day ~ 1 + variety_type_NMWV + s(harvest_day_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  revenue_wheat ~ 1 + variety_type_NMWV + s(Res_wheat_sow) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lamdiag1 ~ 1,
  lamdiag2 ~ 1,
  lamdiag3 ~ 1,
  lamdiag4 ~ 1,
  lambda12 ~ 1,
  lambda13 ~ 1,
  lambda14 ~ 1,
  lambda23 ~ 1,
  lambda24 ~ +s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lambda34 ~ 1
)

multivariate_geo_sow_MRF_corr_endo_rev <- bamlss(f_rice_wheat_yield_MRF_corr_endo_rev, type = "modified", family = mvnchol_bamlss(k = 4), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo_rev)

plot(multivariate_geo_sow_MRF_corr_endo_rev)

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
plot(multivariate_geo_sow_MRF_corr_endo_rev, pages = 1, spar = FALSE, rug = TRUE)
dev.off()

# Maximum autocorrelation
plot(multivariate_geo_sow_MRF_corr_endo_rev, which = "max-acf" , spar = FALSE, lag =200)

dev.off()


# Traceplots
#par(mar = c(4, 4, 4, 1))
#plot(multivariate_geo_sow_MRF_corr_endo_rev, which = "samples")
#model = "mu1", term = "(Intercept)"


# 95% Credible Interval for Predictions
p <- predict(multivariate_geo_sow_MRF_corr_endo_rev, type = "parameter", FUN = c95)
p=as.data.frame(p)


## More results -------------------------------
nd <- data.frame("District" = unique(Irrig_Rev_rice_wheat$District))

# Focusing on the cross-equation correlations

## Predict for the structured spatial effects.
p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev <- predict(multivariate_geo_sow_MRF_corr_endo_rev, newdata = nd, term = "s(District,id='mrf1')", intercept = FALSE)

p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt_rev <- as.data.frame(p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev)

#write.csv(p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt_rev, "tables/p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt2_rev.csv")

## And the unstructured spatial effect.
p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev <- predict(multivariate_geo_sow_MRF_corr_endo_rev, newdata = nd, term = "s(District,id='re2')", intercept = FALSE)

# Rice sowing spatial equation

plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu1, id = nd$District, main = expression(mu(rice_sowing)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu1, id = nd$District, main = expression(mu(rice_sowing)), title = "Unstructured spatial effect")

# Rice yield spatial equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu2, id = nd$District, main = expression(mu(rice_revenues)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu2, id = nd$District, main = expression(mu(rice_revenues)), title = "Unstructured spatial effect")

# Wheat sowing equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu3, id = nd$District, main = expression(mu(wheat_sowing)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu3, id = nd$District, main = expression(mu(wheat_sowing)), title = "Unstructured spatial effect")

# Wheat yield spatial equation
plotmap(India_aoi_sp_bnd, x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu4, id = nd$District, main = expression(mu(wheat_revenues)), title = "Structured spatial effect")

plotmap(India_aoi_sp_bnd, x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$mu4, id = nd$District, main = expression(mu(wheat_revenues)), title = "Unstructured spatial effect")

# Focusing on the cross-equation correlations

## MRF smooth effect.
plotmap(India_aoi_sp_bnd,
        x = p_str_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$lambda24, id = nd$District,
        main = expression(lambda(rice, wheat)), title = "Structured spatial effect"
)

## Random effects.
plotmap(India_aoi_sp_bnd,
        x = p_unstr_multivariate_geo_sow_MRF_corr_endo_rice_y_rev$lambda24, id = nd$District, main = expression(lambda(rice, wheat)), title = "Unstructured spatial effect"
)

# Rice sowing equation : Non-linear relationships
# s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15)
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu1", term = "s(gw_2017)")
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu1", term = "s(onset_2017)")
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu1", term = "s(monsoon_onset_dev)")
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu1", term = "s(median_onset_82_15)")
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu1", term = "s(sd_onset_82_15)")


# Rice yield equation
# s(Res_rice_sow) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice)

plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu2", term = "s(Res_rice_sow)")

plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu2", term = "s(g_q5305_irrig_times_rice)")

plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu2", term = "s(nperha_rice)")

plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu2", term = "s(p2o5perha_rice)")

# Wheat sowing equation
# s(harvest_day_rice) + s(gw_2018)
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu3", term = "s(harvest_day_rice)")

# (multivariate_geo_sow_MRF_corr_endo, model = "mu3", term = "s(gw_2018)")

# Wheat yield equation
plot(multivariate_geo_sow_MRF_corr_endo_rev, model = "mu4", term = "s(Res_wheat_sow)")


# Fitted values
multivariate_geo_sow_MRF_corr_endo_fitted_values_rev <- multivariate_geo_sow_MRF_corr_endo_rev$fitted

multivariate_geo_sow_MRF_corr_endo_fitted_values_rev <- as.data.frame(multivariate_geo_sow_MRF_corr_endo_fitted_values_rev)

# Merge the fitted results to the data and export

Irrig_Rev_rice_wheat_Mult_Res_rev <- cbind(Irrig_Rev_rice_wheat, multivariate_geo_sow_MRF_corr_endo_fitted_values_rev)

write.csv(Irrig_Rev_rice_wheat_Mult_Res_rev, "tables/Irrig_Rev_rice_wheat_Mult_Res_rev.csv")

summary(lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res_rev))

summary(lm(revenue_rice ~ revenue_wheat, Irrig_Rev_rice_wheat_Mult_Res_rev))

# 
#library(flexplot)

#flexplot(revenue_rice ~ revenue_wheat, data=Irrig_Rev_rice_wheat_Mult_Res_rev,method = "lm", se = F)



plot((lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res_rev)))

Irrig_Rev_rice_wheat_Mult_Res_sp_rev <- Irrig_Rev_rice_wheat_Mult_Res_rev
coordinates(Irrig_Rev_rice_wheat_Mult_Res_sp_rev) <- c("o_largest_plot_gps_longitude", "o_largest_plot_gps_latitude")

# Map the correlations
library(modelsummary)
mean_na <- function(x) mean(x, na.rm = TRUE)
sd_na <- function(x) SD(x, na.rm = TRUE)


Irrig_Rev_rice_wheat_Mult_Res_dist_rev <- datasummary(Heading("District") * District ~ Heading("N obs") * N + Heading("%") * Percent() + lambda24 * (Mean + SD), data = Irrig_Rev_rice_wheat_Mult_Res_rev, output = "data.frame")

library(dplyr)
Irrig_Rev_rice_wheat_Mult_Res_dist_rev <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist_rev, "Mean_Rice_Wheat_Rho" = "Mean")
Irrig_Rev_rice_wheat_Mult_Res_dist_rev <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist_rev, "SD_Rice_Wheat_Rho" = "SD")

Irrig_Rev_rice_wheat_Mult_Res_dist_rev <- subset(Irrig_Rev_rice_wheat_Mult_Res_dist_rev, Irrig_Rev_rice_wheat_Mult_Res_dist_rev$District != "Purnia")


# Bihar and EUP map
# India district Map

library(geodata)

India <- gadm(country = "IND", level = 2, path = "shp")

plot(India)

India_aoi <- subset(India, India$NAME_1 == "Bihar" | India$NAME_2 %in% c("Ballia", "Chandauli", "Deoria", "Ghazipur", "Kushinagar", "Maharajganj", "Mau", "Siddharth Nagar", "Gorakhpur"))

plot(India_aoi)

plot(India_aoi, add = TRUE)

library(sf)

India_aoi_sf <- st_as_sf(India_aoi)
library(mapview)

mapview(India_aoi_sf)

# Dissolve the district polygons to form new polygon of Bihar and EUP
library(sf)
India_aoi_sf_dis <- st_union(India_aoi_sf)
mapview(India_aoi_sf_dis)

###
library(modelsummary)

India_aoi_sf$District <- India_aoi_sf$NAME_2

India_aoi_sf$District[India_aoi_sf$District == "Purba Champaran"] <- "EastChamparan"
India_aoi_sf$District[India_aoi_sf$District == "Pashchim Champaran"] <- "WestChamparan"
India_aoi_sf$District[India_aoi_sf$District == "Bhojpur"] <- "Arah"
India_aoi_sf$District[India_aoi_sf$District == "Ballia"] <- "Balia"
India_aoi_sf$District[India_aoi_sf$District == "Ghazipur"] <- "Gazipur"
India_aoi_sf$District[India_aoi_sf$District == "Siddharth Nagar"] <- "Siddharthnagar"

rice_wheat_yield_rho_dist_sf_rev <- merge(India_aoi_sf, Irrig_Rev_rice_wheat_Mult_Res_dist_rev, by = "District")

rice_wheat_yield_rho_dist_sf_rev$Mean_Rice_Wheat_Rho <- as.numeric(rice_wheat_yield_rho_dist_sf_rev$Mean_Rice_Wheat_Rho)

library(mapview)
mapview(rice_wheat_yield_rho_dist_sf_rev, zcol = "Mean_Rice_Wheat_Rho", layer.name = "Rice wheat equation correlation")

library(sf)
rice_wheat_yield_rho_dist_sf_sp_rev <- as_Spatial(rice_wheat_yield_rho_dist_sf_rev)
library(tmap)
tmap_mode("plot")
rice_wheat_yield_rho_dist_sf_sp_map_rev <- tm_shape(rice_wheat_yield_rho_dist_sf_sp_rev) +
  tm_polygons(col = "Mean_Rice_Wheat_Rho", title = "Rice wheat revenue \n equation correlation", style = "quantile",palette = "YlGn") +
  tm_layout(legend.outside = TRUE)

rice_wheat_yield_rho_dist_sf_sp_map_rev

tmap_save(rice_wheat_yield_rho_dist_sf_sp_map, "figures/rice_wheat_yield_rho_dist_sf_sp_map .png")










### Factors explaining the strength of the correlation ---------------------------

## Explain yield trade offs-----------
# Multivariate with sowing dates as endogenous
f_rice_wheat_yield_MRF_corr_endo_more <- list(
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
  lambda24 ~ Res_rice_sow+Res_wheat_sow+g_q5305_irrig_times_rice+g_q5305_irrig_times + nperha + p2o5perha+rice_duration_class_long+variety_type_NMWV+s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lambda34 ~ 1
)

multivariate_geo_sow_MRF_corr_endo_more <- bamlss(f_rice_wheat_yield_MRF_corr_endo_more, type = "modified", family = mvnchol_bamlss(k = 4), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo_more)

plot(multivariate_geo_sow_MRF_corr_endo_more)


## Explain economic tradeoffs ---------------

f_rice_wheat_yield_MRF_corr_endo_rev_more <- list(
  sowdate_fmt_rice_day ~ 1 + rice_duration_class_long + s(gw_2017) + s(onset_2017) + s(monsoon_onset_dev) + s(median_onset_82_15) + s(sd_onset_82_15) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  revenue_rice ~ 1 + rice_duration_class_long + s(Res_rice_sow) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  sowdate_fmt_wheat_day ~ 1 + variety_type_NMWV + s(harvest_day_rice) + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  revenue_wheat ~ 1 + variety_type_NMWV + s(Res_wheat_sow) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lamdiag1 ~ 1,
  lamdiag2 ~ 1,
  lamdiag3 ~ 1,
  lamdiag4 ~ 1,
  lambda12 ~ 1,
  lambda13 ~ 1,
  lambda14 ~ 1,
  lambda23 ~ 1,
  lambda24 ~ 1+Res_rice_sow+Res_wheat_sow+g_q5305_irrig_times_rice+g_q5305_irrig_times + nperha + p2o5perha+rice_duration_class_long+variety_type_NMWV+s(District, bs = "mrf", xt = list("penalty" = K)) + s(District, bs = "re"),
  lambda34 ~ 1
)

multivariate_geo_sow_MRF_corr_endo_rev_more <- bamlss(f_rice_wheat_yield_MRF_corr_endo_rev_more, type = "modified", family = mvnchol_bamlss(k = 4), data = Irrig_Rev_rice_wheat)

summary(multivariate_geo_sow_MRF_corr_endo_rev_more)

plot(multivariate_geo_sow_MRF_corr_endo_rev_more)








# Point referenced geoadditive specification --------------------

Irrig_Rev_rice_wheat_sp <- SpatialPointsDataFrame(cbind(Irrig_Rev_rice_wheat$Longitude, Irrig_Rev_rice_wheat$Latitude), data = Irrig_Rev_rice_wheat, proj4string = CRS("+proj=longlat +datum=WGS84"))

### Yield analysis --------------------------------------

### Non structural 
f_rice_wheat_yield_MRF_corr_endo_yld_non_struct_pointgeo<- list(
  b_grain_yield_ton_per_ha_rice ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(Longitude,Latitude),
  l_ton_per_hectare ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(Longitude,Latitude),
  lamdiag1 ~ 1+rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(Longitude,Latitude),
  lamdiag2 ~ 1+variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(Longitude,Latitude),
  lambda12 ~ 1+Res_rice_sow+Res_wheat_sow+g_q5305_irrig_times_rice+g_q5305_irrig_times + nperha + p2o5perha+rice_duration_class_long+variety_type_NMWV+s(sowdate_fmt_rice_day)+s(sowdate_fmt_wheat_day)+s(Longitude,Latitude)
)


multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo <- bamlss(f_rice_wheat_yield_MRF_corr_endo_yld_non_struct_pointgeo, type = "modified", family = mvnchol_bamlss(k = 2), data = Irrig_Rev_rice_wheat_sp)

summary(multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo)

plot(multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo)

multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo_fitted=multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo$fitted.values

multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo_fitted=as.data.frame(multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo_fitted)

Irrig_Rev_rice_wheat_sp$yield_non_struct_lambda12 <- predict(multivariate_geo_sow_MRF_corr_endo_yld_non_struct_pointgeo,model = "lambda12")

library(tmap)
tmap_mode("plot")
tm_shape(Irrig_Rev_rice_wheat_sp) +
  tm_dots(col = "yield_non_struct_lambda12", title = "Rice-wheat yield equation \n correlation parameter", 
          style = "quantile",palette = "-Spectral",
          size = 0.5) +
  tm_layout(legend.outside = TRUE)+
  tm_legend(text.size = 1)  


# Revenue analysis -----------------------------------------------

### Non structural ---------------------------------

f_rice_wheat_yield_MRF_corr_endo_rev_non_struct_pointgeo<- list(
  revenue_rice ~ 1 + rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(Longitude,Latitude),
  revenue_wheat ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(Longitude,Latitude),
  lamdiag1 ~ 1+rice_duration_class_long + s(sowdate_fmt_rice_day) + s(g_q5305_irrig_times_rice) + s(nperha_rice) + s(p2o5perha_rice) + s(Longitude,Latitude),
  lamdiag2 ~ 1+variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(Longitude,Latitude),
  lambda12 ~ 1+g_q5305_irrig_times_rice+g_q5305_irrig_times + nperha + p2o5perha+rice_duration_class_long+variety_type_NMWV+s(sowdate_fmt_rice_day)+s(sowdate_fmt_wheat_day)+s(Longitude,Latitude)
)


multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo <- bamlss(f_rice_wheat_yield_MRF_corr_endo_rev_non_struct_pointgeo, type = "modified", family = mvnchol_bamlss(k = 2), data = Irrig_Rev_rice_wheat_sp)

summary(multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo)

plot(multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo)

multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo_fitted=multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo$fitted.values

multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo_fitted=as.data.frame(multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo_fitted)

Irrig_Rev_rice_wheat_sp$revenue_non_struct_lambda12 <- predict(multivariate_geo_sow_MRF_corr_endo_rev_non_struct_pointgeo,model = "lambda12")

library(tmap)
tmap_mode("plot")
tm_shape(Irrig_Rev_rice_wheat_sp) +
  tm_dots(col = "revenue_non_struct_lambda12", title = "Rice-wheat revenue equation \n correlation parameter", 
          style = "quantile",palette = "-Spectral",
          size = 0.5) +
  tm_layout(legend.outside = TRUE)+
  tm_legend(text.size = 1)  


# Compare the models using DIC

DIC(multivariate_geo_sow_MRF_corr_endo,multivariate_geo_sow_MRF_corr_endo_more)

# library(distreg.vis)
# library(bamlss)
# 
# library(mvnchol)
# if (interactive()) {
#   distreg.vis:: vis()
# }