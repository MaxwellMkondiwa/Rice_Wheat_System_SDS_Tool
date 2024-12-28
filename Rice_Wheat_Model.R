library(sp)
library(mgcv)
library(bamlss)

# Multivariate geoadditive model
# remotes::install_git("https://git.uibk.ac.at/c4031039/mvnchol")
# install_github("https://github.com/meteosimon/mvnchol")
library(mvnchol)

# Remove NAs in the monsoon variables
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$onset_2017)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$monsoon_onset_dev)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$median_onset_82_15)))
Irrig_Rev_rice_wheat <- subset(Irrig_Rev_rice_wheat, !(is.na(Irrig_Rev_rice_wheat$sd_onset_82_15)))

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

summary(multivariate_geo_sow_MRF_corr_endo)

nd <- data.frame("District" = unique(Irrig_Rev_rice_wheat$District))

# Focusing on the cross-equation correlations

## Predict for the structured spatial effects.
p_str_multivariate_geo_sow_MRF_corr_endo_rice_y <- predict(multivariate_geo_sow_MRF_corr_endo, newdata = nd, term = "s(District,id='mrf1')", intercept = FALSE)

p_str_multivariate_geo_sow_MRF_corr_endo_rice_ydt <- as.data.frame(p_str_multivariate_geo_sow_MRF_corr_endo_rice_y)


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


summary(lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res))

summary(lm(b_grain_yield_ton_per_ha_rice ~ l_ton_per_hectare, Irrig_Rev_rice_wheat_Mult_Res))



plot((lm(sowdate_fmt_rice_day ~ mu1, Irrig_Rev_rice_wheat_Mult_Res)))

Irrig_Rev_rice_wheat_Mult_Res_sp <- Irrig_Rev_rice_wheat_Mult_Res
coordinates(Irrig_Rev_rice_wheat_Mult_Res_sp) <- c("o_largest_plot_gps_longitude", "o_largest_plot_gps_latitude")

# Map the correlations
# library(modelsummary)
# Irrig_Rev_rice_wheat_Mult_Res_dist <- datasummary(Heading("District") * District ~ Heading("N obs") * N + Heading("%") * Percent() + lambda24 * (Mean + SD), data = Irrig_Rev_rice_wheat_Mult_Res, output = "data.frame")

# library(dplyr)
# Irrig_Rev_rice_wheat_Mult_Res_dist <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist, "Mean_Rice_Wheat_Rho" = "Mean")
# Irrig_Rev_rice_wheat_Mult_Res_dist <- rename(Irrig_Rev_rice_wheat_Mult_Res_dist, "SD_Rice_Wheat_Rho" = "SD")

# Irrig_Rev_rice_wheat_Mult_Res_dist <- subset(Irrig_Rev_rice_wheat_Mult_Res_dist, Irrig_Rev_rice_wheat_Mult_Res_dist$District != "Purnia")

# rice_wheat_yield_rho_dist_sf <- merge(India_aoi_sf, Irrig_Rev_rice_wheat_Mult_Res_dist, by = "District")

# rice_wheat_yield_rho_dist_sf$Mean_Rice_Wheat_Rho <- as.numeric(rice_wheat_yield_rho_dist_sf$Mean_Rice_Wheat_Rho)
# library(mapview)
# mapview(rice_wheat_yield_rho_dist_sf, zcol = "Mean_Rice_Wheat_Rho", layer.name = "Rice wheat equation correlation")
# library(sf)
# rice_wheat_yield_rho_dist_sf_sp <- as_Spatial(rice_wheat_yield_rho_dist_sf)
# library(tmap)
# tmap_mode("view")
# rice_wheat_yield_rho_dist_sf_sp_map <- tm_shape(rice_wheat_yield_rho_dist_sf_sp) +
#     tm_polygons(col = "Mean_Rice_Wheat_Rho", title = "Rice wheat equation correlation", style = "quantile") +
#     tm_layout(legend.outside = TRUE)

# tmap_save(rice_wheat_yield_rho_dist_sf_sp_map, "figures/rice_wheat_yield_rho_dist_sf_sp_map .png")