## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(h3sdm)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(tidymodels)
library(spatialsample)
library(DALEX)
library(themis)

## -----------------------------------------------------------------------------
cr <- cr_outline_c

## -----------------------------------------------------------------------------
bio <- terra::rast(system.file("extdata", "bioclim_current.tif", package = "h3sdm"))
names(bio) <- gsub(".*bio_", "bio", names(bio))

## -----------------------------------------------------------------------------
h7 <- h3sdm_get_grid(cr, res = 7)

## -----------------------------------------------------------------------------
bio_predictors <- h3sdm_extract_num(bio, h7)
predictors <- h3sdm_predictors(bio_predictors) |>
  dplyr::select(h3_address, bio1, bio12, bio15, geometry)

## -----------------------------------------------------------------------------
ggplot() +
  theme_minimal() +
  geom_sf(data = predictors, aes(fill = bio1)) +
  scale_fill_viridis_c(option = "B")

## -----------------------------------------------------------------------------
pres <- h3sdm_pres("Silverstoneia flotator", cr, res = 7, limit = 10000)

## -----------------------------------------------------------------------------
records <- h3sdm_pa(pres, predictors, n_pseudoabs = 300)

## -----------------------------------------------------------------------------
table(records$presence)

## -----------------------------------------------------------------------------
ggplot() +
  theme_minimal() +
  geom_sf(data = records, aes(fill = presence)) +
  scale_fill_manual(
    name = "Presence/Absence",
    values = c("red", "blue"),
    labels = c("Absence", "Presence")
  )

## -----------------------------------------------------------------------------
dat <- h3sdm_data(records, predictors)

## -----------------------------------------------------------------------------
scv <- h3sdm_spatial_cv(dat, v = 5, repeats = 1)
autoplot(scv) + theme_minimal()

## -----------------------------------------------------------------------------
receta <- h3sdm_recipe(dat) |>
  themis::step_downsample(presence)

modelo <- parsnip::logistic_reg() |>
  parsnip::set_engine("glm") |>
  parsnip::set_mode("classification")

## -----------------------------------------------------------------------------
wf <- h3sdm_workflow(modelo, receta)

## -----------------------------------------------------------------------------
presence_data <- dat |>
  dplyr::filter(presence == 1)

f <- h3sdm_fit_model(wf, scv, presence_data)

## -----------------------------------------------------------------------------
evaluation_metrics <- h3sdm_eval_metrics(
  fitted_model  = f$cv_model,
  presence_data = presence_data
)
evaluation_metrics

## -----------------------------------------------------------------------------
ggplot(evaluation_metrics, aes(.metric, mean)) +
  theme_minimal() +
  geom_col(width = 0.03, color = "dodgerblue3", fill = "dodgerblue3") +
  geom_point(size = 3, color = "orange") +
  ylim(0, 1)

## -----------------------------------------------------------------------------
p <- h3sdm_predict(f, predictors)

## -----------------------------------------------------------------------------
ggplot() +
  theme_minimal() +
  geom_sf(data = p, aes(fill = prediction)) +
  scale_fill_viridis_c(name = "Habitat\nsuitability", option = "B", direction = -1)

## -----------------------------------------------------------------------------
e <- h3sdm_explain(f$final_model, data = dat)

predictors_to_evaluate <- setdiff(names(e$data), c("h3_address", "x", "y", "presence"))

var_imp <- DALEX::model_parts(
  explainer = e,
  variables = predictors_to_evaluate
)
plot(var_imp)

## -----------------------------------------------------------------------------
pdp_glm <- ingredients::partial_dependence(e, variables = c("bio12", "bio1", "bio15"))
plot(pdp_glm)

## -----------------------------------------------------------------------------
bio_future <- terra::rast(system.file("extdata", "bioclim_future.tif", package = "h3sdm"))
names(bio_future) <- c("bio1", "bio12", "bio15")

bio_future_predictors <- h3sdm_extract_num(bio_future, h7)
predictors_future <- h3sdm_predictors(bio_future_predictors)

p_future <- h3sdm_predict(f, predictors_future)

## -----------------------------------------------------------------------------
ggplot() +
  theme_minimal() +
  geom_sf(data = p_future, aes(fill = prediction)) +
  scale_fill_viridis_c(name = "Habitat\nsuitability", option = "B", direction = -1)

