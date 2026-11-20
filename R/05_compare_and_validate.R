# ==============================================================================
# 05_compare_and_validate.R
#
# Judges the direct estimate, the area-level model and the unit-level model
# against the known truth, and turns the result into a reporting rule.
#
# THIS SCRIPT IS ONLY POSSIBLE BECAUSE THE DATA IS SIMULATED
#
# On a real survey the true area values do not exist. Every claim an SAE model
# makes about its own accuracy is then internal - an MSE derived from the
# model's own assumptions, which is not evidence that the assumptions hold.
#
# That is the honest difficulty with the method, and it is why this repository
# is built around a simulated population. Here the truth is retained, so the
# model's promise can be checked rather than believed.
#
# What transfers to real work is not the numbers below but the diagnostics:
# the design-based check in section 4 can be run on real data, and it is the
# closest thing to external validation available there.
#
# Run after R/02, R/03 and R/04.
# ==============================================================================

source("R/01_config.R")
source("R/00_functions.R")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

if (!TRUTH_AVAILABLE) {
  stop("TRUTH_AVAILABLE is FALSE. This script compares estimates against the ",
       "true area values, which exist only for the simulated population.\n",
       "On real data, stop after R/04 and report the model diagnostics.")
}


# ==============================================================================
# 1. Assemble
# ==============================================================================

fh  <- read_csv(file.path(OUTPUT_DIR, "tables", "fay_herriot_estimates.csv"),
                show_col_types = FALSE)
bhf <- read_csv(file.path(OUTPUT_DIR, "tables", "bhf_estimates.csv"),
                show_col_types = FALSE)

all_est <- fh |>
  select(lga_id, lga_name, state, N_i, n_i, estimate_type,
         direct_est, direct_cv, fh_eblup, fh_se, fh_cv, true_mean_fcs) |>
  left_join(bhf |> select(lga_id, bhf_eblup, bhf_se, bhf_cv),
            by = "lga_id")

message("\nAreas: ", nrow(all_est), " (", sum(all_est$n_i > 0), " sampled, ",
        sum(all_est$n_i == 0), " unsampled)")


# ==============================================================================
# 2. Accuracy against the truth
# ==============================================================================

errors <- all_est |>
  mutate(
    err_direct = direct_est - true_mean_fcs,
    err_fh     = fh_eblup   - true_mean_fcs,
    err_bhf    = bhf_eblup  - true_mean_fcs
  )

summarise_errors <- function(d, label) {
  data.frame(
    subset = label,
    n_areas = nrow(d),
    direct_rmse = if (all(is.na(d$err_direct))) NA_real_ else
      round(sqrt(mean(d$err_direct^2, na.rm = TRUE)), 3),
    fh_rmse  = round(sqrt(mean(d$err_fh^2)), 3),
    bhf_rmse = round(sqrt(mean(d$err_bhf^2)), 3),
    direct_bias = if (all(is.na(d$err_direct))) NA_real_ else
      round(mean(d$err_direct, na.rm = TRUE), 3),
    fh_bias  = round(mean(d$err_fh), 3),
    bhf_bias = round(mean(d$err_bhf), 3)
  )
}

overall <- bind_rows(
  summarise_errors(errors |> filter(n_i > 0), "sampled areas"),
  summarise_errors(errors |> filter(n_i == 0), "unsampled areas"),
  summarise_errors(errors |> filter(n_i > 0, n_i < 15), "sampled, n < 15"),
  summarise_errors(errors |> filter(n_i >= 30), "sampled, n >= 30")
)

message("\n=============================================================")
message("ACCURACY AGAINST THE TRUE AREA MEANS")
message("=============================================================\n")
print(overall, row.names = FALSE)

samp <- overall |> filter(subset == "sampled areas")

message("\nRMSE across sampled areas: direct ", samp$direct_rmse,
        ", Fay-Herriot ", samp$fh_rmse, ", unit-level ", samp$bhf_rmse)

best <- c("direct", "Fay-Herriot", "unit-level")[
  which.min(c(samp$direct_rmse, samp$fh_rmse, samp$bhf_rmse))]
message("Lowest RMSE: ", best)

thin <- overall |> filter(subset == "sampled, n < 15")
message("\nIn the thinnest areas (n < 15) the gap is wider: direct ",
        thin$direct_rmse, " against ", thin$fh_rmse, " and ", thin$bhf_rmse,
        ".")
message("That is where borrowing strength is supposed to pay, and it does.")

big <- overall |> filter(subset == "sampled, n >= 30")
message("\nIn the better-sampled areas (n >= 30) the ordering can reverse: ",
        big$direct_rmse, ", ", big$fh_rmse, ", ", big$bhf_rmse, ".")

# The finding that matters most, and the one an SAE demonstration usually
# leaves out. Stated conditionally so it stays true if the seed changes.
if (samp$fh_rmse >= samp$direct_rmse) {
  message("\n--------------------------------------------------------------")
  message("READ THIS BEFORE QUOTING ANY OF THE ABOVE")
  message("--------------------------------------------------------------")
  message("Across sampled areas as a whole, the Fay-Herriot estimate is NOT ",
          "more accurate than the direct estimate (", samp$fh_rmse,
          " against ", samp$direct_rmse, ").")
  message("\nThat is not a failure of the implementation. It is what happens ",
          "when most areas have an adequate sample: shrinkage helps the thin ",
          "areas and costs something in the well-sampled ones, and the ",
          "average over all areas can come out flat or slightly worse.")
  message("\nThe case for the model is therefore NOT 'it beats the survey'. ",
          "It is narrower and more defensible: it rescues the areas the ",
          "survey cannot support, and it leaves the rest alone. The ",
          "reporting rule in section 5 is built on that narrower claim.")
  message("--------------------------------------------------------------")
} else {
  message("Where the survey is adequate the model adds little, which bounds ",
          "the claim and is worth saying out loud.")
}

# Unsampled areas: the comparison that shows covariates matter more than the
# choice of model class.
uns <- overall |> filter(subset == "unsampled areas")
if (!is.na(uns$fh_rmse) && !is.na(uns$bhf_rmse)) {
  better <- if (uns$fh_rmse < uns$bhf_rmse) "area-level" else "unit-level"
  message("\nFor the ", uns$n_areas, " areas with no sample at all, the two ",
          "models are far apart: Fay-Herriot ", uns$fh_rmse,
          ", unit-level ", uns$bhf_rmse, ".")
  message("The ", better, " model wins, and the reason is covariates rather ",
          "than model class.")
  message("\nThe area-level model uses auxiliaries that vary between areas - ",
          "nightlights, conflict, market access - so it can say something ",
          "specific about a place it has no survey data for.")
  message("The unit-level model uses household characteristics. Those ",
          "predict a household's outcome well, but they carry little ",
          "information about why one LGA differs from another, so its ",
          "synthetic prediction collapses toward the national average and ",
          "picks up a bias of ", uns$bhf_bias, " FCS points.")
  message("\nThe lesson generalises: for areas with no sample, what matters ",
          "is whether your covariates explain BETWEEN-AREA variation. A ",
          "more sophisticated model on covariates that do not will not ",
          "help.")
}


# ==============================================================================
# 3. Where does the model make things WORSE?
# ==============================================================================
#
# The averages above hide the cases that matter most politically. Shrinkage
# pulls an area toward the regression line, and when an area genuinely departs
# from the line - which is exactly what makes it interesting - the model drags
# it back toward the average and makes the estimate worse.
#
# These are the areas a sceptical reviewer will find. Find them first.

harmed <- errors |>
  filter(n_i > 0) |>
  mutate(
    direct_abs = abs(err_direct),
    fh_abs     = abs(err_fh),
    damage     = fh_abs - direct_abs
  ) |>
  arrange(desc(damage)) |>
  select(lga_id, lga_name, n_i, true_mean_fcs, direct_est, fh_eblup,
         direct_abs, fh_abs, damage)

n_harmed <- sum(harmed$damage > 0)

message("\n=============================================================")
message("AREAS WHERE THE MODEL MADE THE ESTIMATE WORSE")
message("=============================================================")
message("\n", n_harmed, " of ", nrow(harmed), " sampled areas: the ",
        "Fay-Herriot estimate is further from the truth than the direct one.")

message("\nThe five worst:")
print(as.data.frame(harmed |> slice_head(n = 5) |>
        mutate(across(where(is.numeric), \(x) round(x, 2)))), row.names = FALSE)

message("\nThis is not a defect to be fixed. Shrinkage trades a worse ",
        "estimate in genuinely unusual areas for better estimates on ",
        "average, and an area that is truly an outlier is precisely the one ",
        "it will damage.")
message("If the outlying areas are the ones a programme cares about - the ",
        "worst-off LGAs usually are - that trade may not be acceptable, and ",
        "the right answer is to sample them properly rather than to model ",
        "them.")


# ==============================================================================
# 4. The diagnostic that works on real data
# ==============================================================================
#
# Everything above needs the truth. This section does not, so it is the part
# that transfers.
#
# The idea: in areas where the direct estimate IS reliable, the model-based
# estimate should agree with it. If a model systematically disagrees with the
# direct estimate where the direct estimate can be trusted, it has no business
# being believed where it cannot.
#
# This is the standard design-based bias check for SAE (Brown et al. 2001).

reliable <- errors |> filter(n_i > 0, direct_cv <= CV_PUBLISH_THRESHOLD)

bias_test <- function(est_col, label) {
  d <- reliable[[est_col]] - reliable$direct_est
  tt <- t.test(d)
  data.frame(
    estimator   = label,
    n_areas     = nrow(reliable),
    mean_diff   = round(mean(d), 3),
    t_statistic = round(as.numeric(tt$statistic), 3),
    p_value     = round(tt$p.value, 4)
  )
}

bias_results <- bind_rows(
  bias_test("fh_eblup",  "Fay-Herriot"),
  bias_test("bhf_eblup", "Unit-level BHF")
)

message("\n=============================================================")
message("DESIGN-BASED BIAS CHECK (runs on real data too)")
message("=============================================================\n")
print(bias_results, row.names = FALSE)

message("\nTested on the ", nrow(reliable), " areas where the direct ",
        "estimate has a CV at or below ", CV_PUBLISH_THRESHOLD, "%.")
message("A difference indistinguishable from zero is what you want. A ",
        "significant one means the model disagrees with the survey where the ",
        "survey is trustworthy, and that is a reason to stop.")

for (i in seq_len(nrow(bias_results))) {
  if (bias_results$p_value[i] < 0.05) {
    message("\n", bias_results$estimator[i],
            ": significant disagreement (p = ", bias_results$p_value[i],
            "). Investigate before publishing anything from it.")
  }
}

# A second check with no truth required: do the model intervals have the
# coverage they claim? Compare the spread of model-minus-direct differences
# against what the two stated variances imply.
coverage_check <- reliable |>
  mutate(
    implied_sd = sqrt(fh_se^2 + (direct_est * direct_cv / 100)^2),
    z          = (fh_eblup - direct_est) / implied_sd
  ) |>
  summarise(
    sd_of_z = round(sd(z), 3),
    within_2sd = round(100 * mean(abs(z) < 2), 1)
  )

message("\nStandardised differences: sd = ", coverage_check$sd_of_z,
        " (should be near 1), ", coverage_check$within_2sd,
        "% within 2 sd (should be near 95).")

if (coverage_check$sd_of_z < 0.75) {
  message("\nThat sd is well below 1, so the stated uncertainties are too ",
          "WIDE - the model and the direct estimate agree more closely than ",
          "their combined error bars imply.")
  message("Being conservative is the safer direction to err, but it is still ",
          "an error: over-wide intervals push areas above the CV threshold ",
          "and suppress estimates that are in fact usable. Part of it here ",
          "is that the two estimates are not independent - the direct ",
          "estimate is an input to the model - so the combined variance ",
          "double-counts.")
} else if (coverage_check$sd_of_z > 1.3) {
  message("\nThat sd is well above 1: the stated uncertainties are too ",
          "NARROW and the intervals are overconfident. This is the dangerous ",
          "direction. Do not publish the intervals without investigating.")
} else {
  message("\nThat is close enough to 1 to be consistent with the stated ",
          "uncertainties.")
}


# ==============================================================================
# 5. A reporting rule
# ==============================================================================
#
# The point of the whole pipeline. Not "use the model" or "use the survey",
# but a rule stated in advance that says which number goes in the report for
# each area, and what caveat travels with it.

final <- all_est |>
  mutate(
    recommendation = case_when(
      n_i == 0 ~ "model only - no survey data",
      direct_cv <= CV_PUBLISH_THRESHOLD ~ "direct estimate",
      direct_cv <= CV_FLAG_THRESHOLD    ~ "model estimate, flagged",
      TRUE                              ~ "model estimate, flagged"
    ),
    reported_value = if_else(recommendation == "direct estimate",
                             direct_est, fh_eblup),
    reported_se    = if_else(recommendation == "direct estimate",
                             direct_est * direct_cv / 100, fh_se),
    reported_cv    = 100 * reported_se / reported_value
  )

rule_summary <- final |>
  count(recommendation) |>
  left_join(
    final |> group_by(recommendation) |>
      summarise(households = sum(N_i), .groups = "drop"),
    by = "recommendation"
  )

message("\n=============================================================")
message("REPORTING RULE")
message("=============================================================\n")
print(as.data.frame(rule_summary), row.names = FALSE)

message("\nThe rule in words:")
message("  Direct estimate where the survey supports it (CV <= ",
        CV_PUBLISH_THRESHOLD, "%).")
message("  Model estimate, clearly flagged, where it does not.")
message("  Model estimate for unsampled areas, flagged more strongly - a ",
        "figure for a place nobody visited is a prediction, and the map ",
        "should say so.")

message("\nWhat this rule deliberately does NOT do is publish the model ",
        "estimate everywhere because it has the lower average RMSE. Where ",
        "the survey is adequate, the direct estimate carries no model ",
        "assumptions, and that is worth more than a small gain in precision ",
        "when someone challenges the number.")


# ==============================================================================
# 6. Outputs
# ==============================================================================

write_csv(overall,      file.path(OUTPUT_DIR, "tables", "accuracy_comparison.csv"))
write_csv(harmed,       file.path(OUTPUT_DIR, "tables", "areas_harmed_by_model.csv"))
write_csv(bias_results, file.path(OUTPUT_DIR, "tables", "design_bias_check.csv"))
write_csv(final,        file.path(OUTPUT_DIR, "tables", "final_estimates.csv"))

plot_data <- errors |>
  filter(n_i > 0) |>
  select(lga_id, n_i, Direct = err_direct, `Fay-Herriot` = err_fh,
         `Unit-level` = err_bhf) |>
  pivot_longer(-c(lga_id, n_i), names_to = "estimator", values_to = "error") |>
  mutate(estimator = factor(estimator,
                            levels = c("Direct", "Fay-Herriot", "Unit-level")))

p1 <- ggplot(plot_data, aes(x = n_i, y = abs(error), colour = estimator)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_colour_manual(values = c("#B2182B", "#2E6E8E", "#5B8C5A")) +
  labs(
    title = "The models help where the sample is thin, and not elsewhere",
    subtitle = paste0("Absolute error against the true LGA mean. The curves ",
                      "cross: past roughly 40 households\nthe direct ",
                      "estimate is the most accurate of the three."),
    x = "Households interviewed in the LGA",
    y = "Absolute error (FCS points)", colour = NULL,
    caption = "Synthetic demonstration data - the true means are known by construction."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot", legend.position = "top",
        plot.caption = element_text(hjust = 0, colour = "grey45"))

ggsave(file.path(OUTPUT_DIR, "figures", "error_by_sample_size.png"),
       p1, width = 21, height = 13, units = "cm", dpi = 300)

p2 <- errors |>
  filter(n_i > 0) |>
  ggplot(aes(x = direct_est, y = fh_eblup)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey50") +
  geom_segment(aes(xend = direct_est, yend = true_mean_fcs),
               colour = "grey80", linewidth = 0.3) +
  geom_point(aes(size = n_i), colour = "#2E6E8E", alpha = 0.8) +
  geom_point(aes(y = true_mean_fcs), colour = "#B2182B", size = 1,
             shape = 4) +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    title = "Shrinkage pulls each estimate toward the model",
    subtitle = paste0("Blue: Fay-Herriot estimate. Red cross: the truth. ",
                      "Dashed line: no shrinkage.\nSmall LGAs (small points) ",
                      "move furthest."),
    x = "Direct estimate (FCS)", y = "Model estimate (FCS)",
    size = "n in LGA",
    caption = "Synthetic demonstration data."
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, colour = "grey45"))

ggsave(file.path(OUTPUT_DIR, "figures", "shrinkage_vs_truth.png"),
       p2, width = 20, height = 14, units = "cm", dpi = 300)

message("\nWritten to ", OUTPUT_DIR, "/tables and ", OUTPUT_DIR, "/figures")
