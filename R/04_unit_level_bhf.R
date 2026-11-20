# ==============================================================================
# 04_unit_level_bhf.R
#
# The unit-level model. Works on individual household records rather than on
# area summaries, and is more precise than the area-level model when its
# requirements are met.
#
# THE REQUIREMENT, STATED FIRST BECAUSE IT DECIDES EVERYTHING
#
# The Battese-Harter-Fuller predictor needs the POPULATION mean of every
# covariate in every area:
#
#   theta_i = Xbar_i'beta + gamma_i (ybar_i - xbar_i'beta)
#
# Xbar_i cannot come from the survey. If the survey could supply a reliable
# population mean for the area, the area would not be small and none of this
# would be necessary. It has to come from a census, a registry, or an
# administrative dataset.
#
# And the covariates must mean the same thing in both sources. A household
# roster counted one way in your questionnaire and another way in the census is
# not the same variable; using it as though it were introduces a bias that no
# diagnostic in this script detects.
#
# That requirement, not the mathematics, is why area-level models are far more
# common in humanitarian work. Census microdata in the settings this pipeline
# is written for is often a decade old, contested, or absent.
#
# Run after R/02.
# ==============================================================================

source("R/01_config.R")
source("R/00_functions.R")

library(dplyr)
library(readr)
library(ggplot2)

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required.\ninstall.packages('lme4')")
}


# ==============================================================================
# 1. Load
# ==============================================================================

sample_df <- read_csv(file.path(DATA_DIR, "survey_sample.csv"),
                      show_col_types = FALSE)
pop_means <- read_csv(file.path(DATA_DIR, "population_covariate_means.csv"),
                      show_col_types = FALSE)
areas     <- read_csv(file.path(DATA_DIR, "area_frame.csv"),
                      show_col_types = FALSE)

pop_sizes <- areas |> select(lga_id, N_i)

message("\nSample           : ", nrow(sample_df), " households")
message("Population means : ", nrow(pop_means), " areas")
message("Covariates       : ", paste(BHF_COVARIATES, collapse = ", "))


# ==============================================================================
# 2. Check the two sources agree
# ==============================================================================
#
# Before fitting anything, compare the sample mean of each covariate against
# the population mean from the register. They should differ by sampling error
# and nothing else.
#
# A large, consistent gap means one of two things, and both are serious: the
# sample is not representative, or the variable is not measured the same way in
# the two sources. Either invalidates the model, and neither shows up in the
# fit statistics.

comparability <- sample_df |>
  group_by(lga_id) |>
  summarise(across(all_of(BHF_COVARIATES), mean, .names = "samp_{.col}"),
            n_i = n(), .groups = "drop") |>
  left_join(pop_means, by = "lga_id")

gaps <- sapply(BHF_COVARIATES, function(v) {
  mean(comparability[[paste0("samp_", v)]] - comparability[[v]], na.rm = TRUE)
})

rel_gaps <- sapply(BHF_COVARIATES, function(v) {
  100 * mean((comparability[[paste0("samp_", v)]] - comparability[[v]]) /
               comparability[[v]], na.rm = TRUE)
})

check <- data.frame(
  covariate     = BHF_COVARIATES,
  mean_gap      = round(as.numeric(gaps), 4),
  pct_gap       = round(as.numeric(rel_gaps), 2)
)

message("\nSample mean minus population mean, averaged over areas:")
print(check, row.names = FALSE)

flagged <- check |> filter(abs(pct_gap) > 10)
if (nrow(flagged) > 0) {
  message("\n", nrow(flagged), " covariate(s) differ by more than 10% between ",
          "the survey and the population source. Investigate before using ",
          "this model - a definitional mismatch here biases every estimate ",
          "and will not appear in any fit statistic.")
} else {
  message("\nNo covariate differs by more than 10%. The two sources appear ",
          "comparable, which is a necessary condition for this model and not ",
          "a sufficient one.")
}


# ==============================================================================
# 3. Fit
# ==============================================================================

bhf_formula <- reformulate(BHF_COVARIATES, response = "fcs")

fit <- lme4::lmer(
  stats::update(bhf_formula, . ~ . + (1 | lga_id)),
  data = sample_df, REML = TRUE
)

vc <- as.data.frame(lme4::VarCorr(fit))
sigma_v2 <- vc$vcov[vc$grp == "lga_id"]
sigma_e2 <- vc$vcov[vc$grp == "Residual"]

message("\n=============================================================")
message("UNIT-LEVEL MODEL FIT")
message("=============================================================")
message("Between-area variance  (sigma_v2): ", round(sigma_v2, 3))
message("Within-area variance   (sigma_e2): ", round(sigma_e2, 3))
message("Intra-area correlation          : ",
        round(sigma_v2 / (sigma_v2 + sigma_e2), 4))

message("\nFixed effects:")
print(round(lme4::fixef(fit), 4))

message("\nThe intra-area correlation is what makes area estimation possible ",
        "at all. If it were zero, households in an area would be no more ",
        "alike than households anywhere, and there would be nothing ",
        "area-specific to predict.")


# ==============================================================================
# 4. Predict
# ==============================================================================

bhf <- fit_bhf(
  formula   = bhf_formula,
  data      = sample_df,
  area_var  = "lga_id",
  pop_means = pop_means,
  pop_sizes = pop_sizes
)

sampled_bhf <- bhf |> filter(!is.na(n_i))

message("\nShrinkage weights (gamma): min ", round(min(sampled_bhf$gamma), 3),
        ", median ", round(median(sampled_bhf$gamma), 3),
        ", max ", round(max(sampled_bhf$gamma), 3))

message("\nThese gammas are higher than the area-level model's, and the ",
        "reason is structural: the unit-level model uses every household ",
        "record rather than one summary per area, so there is more ",
        "information behind each area's own contribution.")


# ==============================================================================
# 5. Unsampled areas
# ==============================================================================
#
# With no sample in an area, ybar_i does not exist and the predictor reduces to
# the synthetic part, Xbar_i'beta. Note what this depends on: the population
# covariate means for an area nobody visited. If the register covers it, the
# estimate is available; if the register is as incomplete as the survey, the
# unit-level model offers nothing the area-level model does not.

unsampled_bhf <- bhf |> filter(is.na(n_i))

if (nrow(unsampled_bhf) > 0) {

  # fit_bhf() already computed Xbar_i'beta for every area in pop_means, which
  # is exactly the synthetic predictor. For an unsampled area that is the whole
  # estimate - there is no ybar_i to correct it with.
  #
  # The variance of a purely synthetic prediction includes sigma_v2, because
  # the area's own departure is unknown. Omitting it is the standard way these
  # estimates get published looking more precise than they are.
  unsampled_bhf <- unsampled_bhf |>
    mutate(
      eblup = synthetic,
      gamma = 0,
      mse   = sigma_v2
    )

  message("\n", nrow(unsampled_bhf), " unsampled LGAs predicted from the ",
          "population covariate means alone.")
  message("Their MSE is sigma_v2 itself - the full between-area variance, ",
          "because nothing is known about their individual departure.")

  bhf <- bind_rows(sampled_bhf, unsampled_bhf)
} else {
  bhf <- sampled_bhf
}

bhf <- bhf |>
  mutate(
    bhf_se = sqrt(mse),
    bhf_cv = cv_pct(eblup, mse),
    estimate_type = if_else(is.na(n_i) | n_i == 0,
                            "synthetic (no sample)", "model (sampled)")
  )


# ==============================================================================
# 6. A caveat about the outcome type
# ==============================================================================
#
# This model is linear and assumes normal errors. The food consumption score is
# continuous and roughly symmetric, so that is defensible.
#
# It would NOT be defensible for a binary indicator - the proportion of
# households with poor food consumption, say - which is what humanitarian
# reporting more often asks for. A linear model on a 0/1 outcome produces
# predictions outside [0,1] and variance estimates that are wrong by
# construction.
#
# The correct tool there is a generalised linear mixed model, and the
# small-area literature on it (Molina & Rao's empirical best prediction for
# poverty indicators) is a larger subject than this repository covers. The
# honest position is to say so rather than to run the linear model on a
# proportion and hope.

message("\nNOTE: this model is linear. It suits the continuous FCS and would ",
        "NOT suit a binary indicator such as the proportion below a ",
        "threshold - see docs/methodology.md.")


# ==============================================================================
# 7. Outputs
# ==============================================================================

out <- bhf |>
  left_join(areas |> select(lga_id, lga_name, state, any_of("true_mean_fcs")),
            by = "lga_id") |>
  mutate(n_i = coalesce(n_i, 0L)) |>
  select(lga_id, lga_name, state, N_i, n_i, estimate_type,
         sample_mean = ybar, synthetic, gamma,
         bhf_eblup = eblup, bhf_se, bhf_cv,
         any_of("true_mean_fcs")) |>
  arrange(lga_id)

write_csv(out, file.path(OUTPUT_DIR, "tables", "bhf_estimates.csv"))
saveRDS(list(fit = fit, sigma_v2 = sigma_v2, sigma_e2 = sigma_e2),
        file.path(OUTPUT_DIR, "bhf_model.rds"))

message("\nWritten to ", OUTPUT_DIR, "/tables")
