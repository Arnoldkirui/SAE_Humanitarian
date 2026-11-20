# ==============================================================================
# 03_fay_herriot.R
#
# The area-level model. Takes the direct estimates and their variances from
# R/02, borrows strength from area-level auxiliary data, and produces an EBLUP
# for every LGA - including the ones with no sample.
#
# WHAT THE MODEL IS DOING
#
# For each area it forms a weighted average of two things: the direct estimate,
# which is unbiased but noisy, and a regression prediction from the auxiliary
# data, which is precise but may be biased for that particular area. The weight
# is gamma, and it is not a tuning choice - it falls out of how the variance of
# the direct estimate compares with the variance between areas that the
# auxiliaries cannot explain.
#
#   gamma near 1   the direct estimate is reliable; the model barely moves it
#   gamma near 0   the direct estimate is noise; the regression takes over
#
# An area with no sample has gamma = 0 by construction: there is no direct
# estimate to weight, and the prediction is purely synthetic.
#
# Run after R/02.
# ==============================================================================

source("R/01_config.R")
source("R/00_functions.R")

library(dplyr)
library(readr)
library(ggplot2)


# ==============================================================================
# 1. Assemble the area-level dataset
# ==============================================================================

areas  <- read_csv(file.path(DATA_DIR, "area_frame.csv"), show_col_types = FALSE)
direct <- read_csv(file.path(OUTPUT_DIR, "tables", "direct_estimates.csv"),
                   show_col_types = FALSE)

# Transformations applied before modelling. Skewed covariates entered raw give
# the regression a few high-leverage areas; logging them is routine and should
# be stated rather than silently applied.
areas <- areas |>
  mutate(
    log_nightlights = log(nightlights),
    log_pop_density = log(pop_density),
    sqrt_conflict   = sqrt(conflict_events)
  )

# areas already carries a planned n_i from the sampling design; the realised
# one comes from the sample itself. Drop the planned column so the join does
# not silently produce n_i.x and n_i.y.
fh_data <- areas |>
  select(-any_of("n_i")) |>
  left_join(direct |> select(lga_id, n_i, direct_est, direct_var,
                             direct_se, cv),
            by = "lga_id")

sampled   <- fh_data |> filter(!is.na(direct_est))
unsampled <- fh_data |> filter(is.na(direct_est))

message("\nAreas with a direct estimate: ", nrow(sampled))
message("Areas without             : ", nrow(unsampled))


# ==============================================================================
# 2. Check the auxiliary data before trusting it
# ==============================================================================
#
# Two failures to rule out first, because both produce a model that looks
# excellent and is worthless.
#
# The first is a covariate derived from the survey. An area-level model needs
# auxiliaries available everywhere, and anything computed from the sample is
# missing for exactly the areas that matter most.
#
# The second is a covariate so strongly related to the outcome that it is
# effectively the outcome measured another way. That inflates R-squared,
# collapses sigma_v2 toward zero, and produces EBLUPs that are just the
# regression line.

aux_check <- lm(reformulate(FH_COVARIATES, response = "direct_est"),
                data = sampled)

r2 <- summary(aux_check)$r.squared

message("\nAuxiliary model R-squared on the direct estimates: ", round(r2, 3))

if (r2 > 0.95) {
  message("\nThat is suspiciously high. Check whether one of the covariates ",
          "is derived from the outcome - at this level of fit the model is ",
          "unlikely to be learning anything the auxiliaries did not already ",
          "contain.")
} else if (r2 < 0.2) {
  message("\nThe auxiliaries explain very little. Borrowing strength from ",
          "them will barely move the direct estimates, and for unsampled ",
          "areas the synthetic prediction will be close to the overall mean. ",
          "Consider whether better covariates are available before ",
          "proceeding.")
}

# Multicollinearity. Correlated auxiliaries do not bias the EBLUP but they make
# individual coefficients uninterpretable, and someone will ask what the
# nightlights coefficient means.
cors <- cor(sampled[, FH_COVARIATES], use = "complete.obs")
max_cor <- max(abs(cors[upper.tri(cors)]))
message("Largest pairwise correlation among auxiliaries: ", round(max_cor, 3))
if (max_cor > 0.8) {
  message("Strongly correlated covariates - do not interpret the individual ",
          "coefficients, and expect them to be unstable across samples.")
}


# ==============================================================================
# 3. Fit
# ==============================================================================

X_sampled <- model.matrix(reformulate(FH_COVARIATES), data = sampled)

fh <- fit_fay_herriot(
  y   = sampled$direct_est,
  X   = X_sampled,
  psi = sampled$direct_var
)

message("\n=============================================================")
message("FAY-HERRIOT FIT")
message("=============================================================")
message("sigma_v2 (unexplained between-area variance): ",
        round(fh$sigma_v2, 4))
message("Residual sd between areas                  : ",
        round(sqrt(fh$sigma_v2), 3), " FCS points")

if (fh$sigma_at_zero) {
  message("\nsigma_v2 estimated at zero. Every EBLUP is the regression ",
          "prediction, and the survey is contributing nothing. Investigate ",
          "the auxiliaries before using this.")
}

coefs <- data.frame(
  term     = colnames(X_sampled),
  estimate = round(fh$beta, 4)
)
message("\nCoefficients:")
print(coefs, row.names = FALSE)

message("\nShrinkage weights (gamma): min ", round(min(fh$gamma), 3),
        ", median ", round(median(fh$gamma), 3),
        ", max ", round(max(fh$gamma), 3))
message("A gamma of ", round(median(fh$gamma), 2),
        " at the median means the typical LGA's estimate is about ",
        round(100 * median(fh$gamma)), "% direct estimate and ",
        round(100 * (1 - median(fh$gamma))), "% model.")


# ==============================================================================
# 4. EBLUPs for sampled areas
# ==============================================================================

sampled_out <- sampled |>
  mutate(
    fh_eblup = fh$eblup,
    fh_mse   = fh$mse,
    fh_se    = sqrt(fh$mse),
    fh_cv    = cv_pct(fh_eblup, fh$mse),
    gamma    = fh$gamma,
    synthetic = fh$x_beta,
    shift    = fh_eblup - direct_est,
    estimate_type = "model (sampled)"
  )


# ==============================================================================
# 5. Synthetic estimates for unsampled areas
# ==============================================================================
#
# With no direct estimate, the prediction is x'beta alone. Its variance has two
# parts: the uncertainty in beta, and sigma_v2 itself - because the area's own
# departure v_i is entirely unknown and could be anything the distribution of
# departures allows.
#
# That second part is what makes these estimates much less precise than they
# look, and it is routinely omitted when synthetic estimates are published.

if (nrow(unsampled) > 0) {
  X_unsampled <- model.matrix(reformulate(FH_COVARIATES), data = unsampled)

  XtVinvX_inv <- solve(crossprod(X_sampled,
                                 X_sampled / (fh$sigma_v2 + sampled$direct_var)))

  synth_pred <- as.vector(X_unsampled %*% fh$beta)

  var_beta_part <- sapply(seq_len(nrow(X_unsampled)), function(i) {
    xi <- X_unsampled[i, , drop = FALSE]
    as.numeric(xi %*% XtVinvX_inv %*% t(xi))
  })

  synth_mse <- var_beta_part + fh$sigma_v2

  unsampled_out <- unsampled |>
    mutate(
      fh_eblup = synth_pred,
      fh_mse   = synth_mse,
      fh_se    = sqrt(synth_mse),
      fh_cv    = cv_pct(fh_eblup, synth_mse),
      gamma    = 0,
      synthetic = synth_pred,
      shift    = NA_real_,
      estimate_type = "synthetic (no sample)"
    )

  message("\nSynthetic estimates for ", nrow(unsampled), " unsampled LGAs:")
  message("  CV range: ", round(min(unsampled_out$fh_cv), 1), "% to ",
          round(max(unsampled_out$fh_cv), 1), "%")
  message("  Of that uncertainty, ",
          round(100 * mean(fh$sigma_v2 / synth_mse)),
          "% comes from the unknown area effect alone - the part no amount ",
          "of auxiliary data removes.")

  fh_results <- bind_rows(sampled_out, unsampled_out)
} else {
  fh_results <- sampled_out
}


# ==============================================================================
# 6. Did it help?
# ==============================================================================

comparison <- fh_results |>
  filter(estimate_type == "model (sampled)") |>
  summarise(
    median_direct_cv = round(median(cv), 1),
    median_fh_cv     = round(median(fh_cv), 1),
    n_improved       = sum(fh_cv < cv),
    n_worse          = sum(fh_cv >= cv),
    max_shift        = round(max(abs(shift)), 2)
  )

message("\n=============================================================")
message("PRECISION, DIRECT VERSUS MODEL")
message("=============================================================")
print(as.data.frame(comparison))

message("\nThe MSE is lower than the direct variance in ",
        comparison$n_improved, " of ", nrow(sampled), " sampled LGAs.")
message("Largest movement from the direct estimate: ",
        comparison$max_shift, " FCS points.")
message("\nA reduction in MSE is what the model promises, and it is a ",
        "model-based promise: it holds if the linking model is right. R/05 ",
        "checks it against the truth instead of taking it on trust.")

n_now_publishable <- sum(fh_results$fh_cv <= CV_PUBLISH_THRESHOLD)
message("\nLGAs below the ", CV_PUBLISH_THRESHOLD, "% CV threshold: ",
        n_now_publishable, " of ", nrow(fh_results),
        " (direct estimates managed ", sum(direct$cv <= CV_PUBLISH_THRESHOLD),
        " of ", nrow(areas), ")")


# ==============================================================================
# 7. Outputs
# ==============================================================================

out <- fh_results |>
  mutate(n_i = coalesce(n_i, 0L)) |>
  select(lga_id, lga_name, state, N_i, n_i, estimate_type,
         direct_est, direct_se, direct_cv = cv,
         synthetic, gamma, fh_eblup, fh_se, fh_cv,
         any_of("true_mean_fcs")) |>
  arrange(lga_id)

write_csv(out, file.path(OUTPUT_DIR, "tables", "fay_herriot_estimates.csv"))
saveRDS(fh, file.path(OUTPUT_DIR, "fh_model.rds"))

p <- fh_results |>
  filter(estimate_type == "model (sampled)") |>
  ggplot(aes(x = n_i, y = gamma)) +
  geom_point(colour = "#2E6E8E", size = 2.2, alpha = 0.85) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "How much the model trusts each LGA's own data",
    subtitle = paste0("Shrinkage weight against sample size. At gamma = 1 the ",
                      "estimate is purely the survey;\nat gamma = 0 it is ",
                      "purely the model."),
    x = "Households interviewed in the LGA",
    y = "Shrinkage weight (gamma)",
    caption = "Synthetic demonstration data."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, colour = "grey45"))

ggsave(file.path(OUTPUT_DIR, "figures", "fh_shrinkage.png"),
       p, width = 20, height = 12, units = "cm", dpi = 300)

message("\nWritten to ", OUTPUT_DIR, "/tables")
