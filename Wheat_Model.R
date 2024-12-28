# Wheat model

f_wheat_yield_MRF <- list(
  l_ton_per_hectare ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re"),
  sigma ~ 1 + variety_type_NMWV + s(sowdate_fmt_wheat_day) + g_q5305_irrig_times + nperha + p2o5perha + s(District, bs = "mrf", xt = list("penalty" = K)) +
    s(District, bs = "re")
)
set.seed(321)
b_wheat_yield_MRF <- bamlss(f_wheat_yield_MRF, data = Irrig_Rev_rice_wheat, family = "gaussian")