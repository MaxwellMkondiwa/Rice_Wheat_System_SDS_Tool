###############################################################################
# Interactively visualizing distributional regression models with distreg.vis #
# Paper supplement                                                            #
# File: replicate.R                                                         #
# Description: Replicate the R output and graphs                              #
###############################################################################

## Load all necessary libraries
library("bamlss")
library("gamlss")
library("ggplot2")
library("ISLR")
library("distreg.vis")
library("patchwork")

## Set seed
set.seed(123)

#### ---- MAIN PAPER ---- ####

## --- Section: Motivational Example --- ##

# Data
Wage <- ISLR::Wage
colnames(Wage)[colnames(Wage) == "race"] <- "ethnicity"

# Model Building
wage_model <- bamlss(
  list(
    wage ~ s(age) + ethnicity + year + education,
    sigma ~ s(age) + ethnicity + year + education
  ),
  data = Wage,
  family = lognormal_bamlss()
)

# Model Output
print(summary(wage_model), digits = 1)

# Data.frame for prediction
df <- set_mean(
  model_data(wage_model),
  vary_by = "education"
)
row.names(df) <- levels(Wage$education)
df

# Predict Parameters
pp <- preds(model = wage_model,
            newdata = df)
pp

pp_plot <- plot_dist(wage_model, pred_params = pp)

# Look at influence of age
moments_plot_age <- plot_moments(
  wage_model,
  int_var = "age",
  pred_data = df,
  rug = TRUE,
  samples = TRUE,
  palette = "viridis",
  uncertainty = TRUE
)

# Look at influence on external function
gini <- function(par) {
  2 * pnorm((par[["sigma"]] / 2) * sqrt(2)) - 1
}
moments_plot_exfun <- plot_moments(
  wage_model,
  int_var = "age",
  pred_data = df,
  samples = TRUE,
  uncertainty = TRUE,
  ex_fun = "gini",
  rug = TRUE)

# Look at influence of ethnicity
moments_plot_ethno <- plot_moments(
  wage_model,
  int_var = "ethnicity",
  pred_data = df,
  samples = TRUE,
  palette = "viridis",
  uncertainty = TRUE
)

# Save plot
ggsave(pp_plot, file = "pp_plot.pdf", width = 7.5, height = 4.5, dpi = 320)
ggsave(moments_plot_age, file = "moments_plot_age.pdf", width = 7.5, height = 4.5, dpi = 320)
ggsave(moments_plot_exfun, file = "moments_plot_exfun.pdf", width = 7.5, height = 4.5, dpi = 320)
ggsave(moments_plot_ethno, file = "moments_plot_ethno.pdf", width = 7.5, height = 4.5, dpi = 320)

## --- Section: An Introduction to the Graphical User Interface --- ##

# Opening of App
if (interactive())
  vis()

#### ---- APPENDIX ---- ####

## --- Section: Implementation of distreg.vis --- ##

### Showcase of tools to choose covariate combinations ###
df <- set_mean(input = model_data(wage_model), vary_by = "education")
row.names(df) <- levels(Wage$education)
df

### Showcase of Preds ###
pp_samples <- preds(model = wage_model,
                    newdata = df,
                    what = "samples")
lapply(pp_samples, head, 2)

### Showcase of dists ###
str(dists)

### Showcase of plot_dist() ###

## BE (Beta)
art_data <- model_fam_data(fam_name = "BE")
form_be <- as.formula("BE ~ norm2 + binomial1")
model_be <- gamlss(form_be, sigma.formula = ~ .,
                   data = art_data, family = "BE", trace = FALSE)
ndata_be <- art_data[sample(seq_len(nrow(art_data)), 5),
                     !colnames(art_data) %in% "BE"]
pred_params_be <- preds(model_be, newdata = ndata_be)

## Beta - Plots
plot_dist_beta_pdf <- plot_dist(model_be, pred_params_be, rug = TRUE) +
  labs(title = "b)") + theme(plot.title = element_text(hjust = 0.5)) # pdf
plot_dist_beta_cdf <- plot_dist(model_be, pred_params_be, rug = TRUE, type = "cdf") +
  labs(title = "c)")  + # cdf
  theme(plot.title = element_text(hjust = 0.5))

## GEOM (Geometrical)
art_data <- model_fam_data(fam_name = "GEOM")
form_geom <- as.formula("GEOM ~ norm2 + binomial1")
model_geom <- gamlss(form_geom, sigma.formula = ~ .,
                     data = art_data, family = "GEOM", trace = FALSE)
ndata_geom <- art_data[sample(seq_len(nrow(art_data)), 5),
                       !colnames(art_data) %in% "GEOM"]
pred_params_geom <- preds(model_geom, newdata = ndata_geom)

## GEOM - Plots
plot_dist_geom_pdf <- plot_dist(model_geom, pred_params_geom, rug = TRUE) +
  lims(y = c(NA, 1)) +
  labs(title = "d)") +
  theme(plot.title = element_text(hjust = 0.5)) # pdf
plot_dist_geom_cdf <- plot_dist(model_geom, pred_params_geom, rug = TRUE, type = "cdf") +
  labs(title = "e)") + # cdf
  lims(y = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5))

## Multinomial
art_data_multinom <- model_fam_data(fam_name = "multinomial")
fam_called <- multinomial_bamlss()
model_multinom <- bamlss(list(multinomial ~ norm2 + binomial1,
                              ~ norm2 + binomial1),
                         data = art_data_multinom,
                         family = fam_called, verbose = FALSE)
ndata_multinom <- art_data_multinom[sample(seq_len(nrow(art_data_multinom)), 5),
                                    !colnames(art_data_multinom) %in% "multinomial"]
pred_params_multinom <- preds(model_multinom, newdata = ndata_multinom)

## Multinomial - Plot
plot_dist_multinomial <- plot_dist(model_multinom, pred_params_multinom, rug = TRUE) +
  labs(title = "a)") + # pdf
  theme(plot.title = element_text(hjust = 0.5))

## Pierce all plots together
plot_dist_showcase <- plot_dist_multinomial + {
  plot_dist_beta_pdf +
    plot_dist_beta_cdf +
    plot_dist_geom_pdf +
    plot_dist_geom_cdf
} + plot_layout(ncol = 1, heights = c(1, 3))
ggsave(plot_dist_showcase, file = "plot_dist_showcase.pdf", width = 7.5, height = 7.5, dpi = 320)

### Showcase of plot_moments() ###
df_new <- df[c(1, 3, 5), ]
gini <- function(par) {
  2 * pnorm((par[["sigma"]] / 2) * sqrt(2)) - 1
}
p_1 <- plot_moments(wage_model,
                    int_var = "age",
                    pred_data = df_new,
                    samples = TRUE,
                    uncertainty = TRUE,
                    ex_fun = "gini",
                    rug = TRUE) +
  labs(title = "a)") +
  theme(plot.title = element_text(hjust = 0.5))
p_2 <- plot_moments(wage_model,
                    int_var = "ethnicity",
                    pred_data = df_new,
                    samples = TRUE,
                    uncertainty = TRUE,
                    ex_fun = "gini") +
  labs(title = "b)") +
  theme(plot.title = element_text(hjust = 0.5))

## Pierce plots together
plotmomentsgraph <- p_1 + p_2 + plot_layout(ncol = 1)

### Save both plots
ggsave(plotmomentsgraph, file = "plotmomentsgraph.pdf", width = 10, height = 7.5, dpi = 320)

## --- Section: The Graphical User Interface --- ##

# Opening of App
if (interactive())
  vis()

## --- Section: Additional Graphs --- ##

## Multinomial family
ndata_multinom <- ndata_multinom[c(1, 2, 3), ]
rownames(ndata_multinom) <- c("First Scenario", "Second Scenario", "Third Scenario")
p_multinomial <- plot_moments(model_multinom,
                              int_var = "norm2",
                              pred_data = ndata_multinom)

ggsave(p_multinomial, file = "pmultinomial.pdf", width = 10, height = 4, dpi = 320)
