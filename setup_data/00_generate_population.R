# ==============================================================================
# 00_generate_population.R
#
# Builds a synthetic finite population of households across 60 LGAs, then draws
# a realistic two-stage survey sample from it.
#
# WHY A SIMULATED POPULATION RATHER THAN A SAMPLE DATASET
#
# Every small area estimation method claims to beat the direct estimate. On
# real data that claim cannot be checked, because the true area value is
# exactly what you do not have - if you knew it, you would not be modelling.
#
# Simulating the population first means the truth is retained. R/05 can then
# ask the only question that matters: did the model-based estimate land closer
# to the truth than the direct one, and in which areas did it not?
#
# NOTHING HERE IS REAL. No LGA, household or figure in this dataset
# corresponds to a real one. The structure is shaped like a Nigerian MSNA -
# three states, uneven LGA sample sizes, a food consumption outcome - because
# that is the setting the method is being demonstrated for.
#
# Fixed seed, so every run reproduces exactly.
# ==============================================================================

# ------------------------------------------------------------------------------
# EDIT THIS BLOCK ONLY
# ------------------------------------------------------------------------------

OUT_DIR <- "data"

N_LGAS          <- 60      # areas
TARGET_SAMPLE   <- 2400    # total households sampled
CLUSTER_SIZE    <- 6       # households interviewed per cluster
N_UNSAMPLED_LGAS <- 9      # LGAs the assessment did not reach at all

# ------------------------------------------------------------------------------
# NOTHING BELOW NEEDS EDITING
# ------------------------------------------------------------------------------

set.seed(2026)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---- 1. Area frame -----------------------------------------------------------
#
# Invented LGA names. They are shaped like Nigerian admin-2 units so the
# examples read naturally, but none of these places exist.

prefixes <- c("Gar", "Kun", "Bir", "Dam", "Mai", "Jak", "Nguru", "Bade",
              "Tarm", "Kaga", "Mobb", "Gubi", "Shan", "Rann", "Dikw",
              "Monz", "Askir", "Hawu", "Lassa", "Uba")
suffixes <- c("a", "wa", "ri", "no", "shi", "la")

lga_names <- unique(as.vector(outer(prefixes, suffixes, paste0)))[seq_len(N_LGAS)]

states <- rep(c("Northern", "Central", "Eastern"), length.out = N_LGAS)

areas <- data.frame(
  lga_id   = sprintf("LGA%02d", seq_len(N_LGAS)),
  lga_name = lga_names,
  state    = states,
  stringsAsFactors = FALSE
)

# Population size per LGA, right-skewed as real admin units are
areas$N_i <- round(rlnorm(N_LGAS, log(1800), 0.55))
areas$N_i <- pmin(pmax(areas$N_i, 400), 6000)

# ---- 2. Area-level auxiliary variables ---------------------------------------
#
# These are the covariates an area-level model uses. All of them are things you
# could obtain WITHOUT a household survey, which is the point: auxiliary data
# must be available for every area, including the ones you did not sample well.
#
#   nightlights      VIIRS radiance, a proxy for economic activity
#   pop_density      persons per km2
#   market_distance  travel time to the nearest functioning market, hours
#   conflict_events  incidents in the preceding 12 months
#   ndvi_anomaly     growing-season vegetation departure from the long-term mean
#
# They are generated with realistic correlation structure. Real auxiliary data
# is correlated - nightlights and population density especially - and a model
# fitted on independent covariates behaves better than it should.

state_effect <- c(Northern = -0.55, Central = 0.15, Eastern = 0.40)

latent_development <- rnorm(N_LGAS, state_effect[areas$state], 0.75)

areas$nightlights     <- round(exp(0.9 + 0.75 * latent_development +
                                     rnorm(N_LGAS, 0, 0.35)), 3)
areas$pop_density     <- round(exp(4.1 + 0.55 * latent_development +
                                     rnorm(N_LGAS, 0, 0.45)))
areas$market_distance <- round(pmax(0.2, 3.1 - 0.85 * latent_development +
                                      rnorm(N_LGAS, 0, 0.55)), 2)
areas$conflict_events <- rpois(N_LGAS, lambda = pmax(0.5, 14 -
                                        4.5 * latent_development))
areas$ndvi_anomaly    <- round(rnorm(N_LGAS, 0.02 + 0.05 * latent_development,
                                     0.12), 3)

# ---- 3. True area means ------------------------------------------------------
#
# The outcome is the Food Consumption Score. The true area mean is driven by
# the auxiliaries PLUS an area-specific departure the auxiliaries cannot
# explain.
#
# That departure is the whole reason SAE is a model rather than a regression.
# If the auxiliaries explained everything, you would predict directly from them
# and never look at the survey. If they explained nothing, borrowing strength
# would add nothing. Real settings sit in between, and sigma_v2 is what the
# model estimates to locate them.

areas$true_fcs_mean <-
  34.5 +
  3.2  * scale(log(areas$nightlights))[, 1] +
  1.6  * scale(log(areas$pop_density))[, 1] +
 -2.4  * scale(areas$market_distance)[, 1] +
 -1.9  * scale(sqrt(areas$conflict_events))[, 1] +
  1.3  * scale(areas$ndvi_anomaly)[, 1] +
  rnorm(N_LGAS, 0, 2.8)                      # the unexplained area departure

areas$true_fcs_mean <- round(pmax(10, areas$true_fcs_mean), 3)

# ---- 4. Household population -------------------------------------------------
#
# Household-level covariates are generated for the whole population, because
# the unit-level model needs their POPULATION means - which in real work come
# from a census or register, not from the survey.

make_households <- function(i) {
  n <- areas$N_i[i]

  hh_size <- pmax(1, round(rnorm(n, 6.2 + 0.4 * scale(areas$pop_density)[i], 2.4)))
  displaced <- rbinom(n, 1, plogis(-0.9 + 0.08 * areas$conflict_events[i]))
  female_headed <- rbinom(n, 1, 0.22)
  dependency <- round(pmin(0.95, pmax(0.05, rnorm(n, 0.52, 0.16))), 3)

  fcs <- areas$true_fcs_mean[i] +
    -0.55 * (hh_size - 6.2) +
    -4.10 * displaced +
    -2.30 * female_headed +
    -6.50 * (dependency - 0.52) +
    rnorm(n, 0, 11.5)

  data.frame(
    lga_id        = areas$lga_id[i],
    hh_id         = sprintf("%s_H%05d", areas$lga_id[i], seq_len(n)),
    hh_size       = hh_size,
    displaced     = displaced,
    female_headed = female_headed,
    dependency    = dependency,
    fcs           = round(pmax(0, pmin(112, fcs)), 2),
    stringsAsFactors = FALSE
  )
}

population <- do.call(rbind, lapply(seq_len(N_LGAS), make_households))

# Recompute the true area mean from the realised population. This is the
# finite-population truth the estimates are judged against - slightly different
# from the model mean above, which is correct: a finite population is itself a
# realisation.
truth <- aggregate(fcs ~ lga_id, population, mean)
names(truth)[2] <- "true_mean_fcs"
truth$true_mean_fcs <- round(truth$true_mean_fcs, 3)

areas <- merge(areas, truth, by = "lga_id")

# ---- 5. Draw the sample ------------------------------------------------------
#
# Two-stage, and deliberately uneven. A national MSNA is powered for STATE
# level; LGA sample sizes then fall out of the cluster allocation and vary by
# an order of magnitude. Some LGAs get 6 households, some get 90.
#
# That unevenness is the entire problem SAE addresses, so it is built in rather
# than smoothed away.

alloc_weight <- sqrt(areas$N_i) * rlnorm(N_LGAS, 0, 0.78)
n_clusters <- pmax(1, round(TARGET_SAMPLE / CLUSTER_SIZE *
                              alloc_weight / sum(alloc_weight)))

# Cap the largest allocation. Real cluster allocations are bounded by what a
# field team can cover in an LGA, and an uncapped draw produces one implausible
# giant that would flatter the direct estimator in that area.
n_clusters <- pmin(n_clusters, 12)

# Some LGAs are not sampled at all. This is the normal situation, not an edge
# case: an assessment powered for state level reaches the LGAs its clusters
# happen to fall in, and access constraints remove others entirely.
#
# These areas are the sharpest test of the whole exercise. They have no direct
# estimate - not an imprecise one, none at all - so the only thing that can say
# anything about them is a model. Whether that is acceptable is a judgement the
# analyst has to make explicitly, and R/05 shows what it costs.
unsampled <- sample(seq_len(N_LGAS), N_UNSAMPLED_LGAS)
n_clusters[unsampled] <- 0
areas$n_clusters <- n_clusters
areas$n_i <- n_clusters * CLUSTER_SIZE

sample_lga <- function(i) {
  if (areas$n_clusters[i] == 0) return(NULL)

  pop_i <- population[population$lga_id == areas$lga_id[i], ]
  n_take <- min(areas$n_i[i], nrow(pop_i))

  # Cluster sampling: households within a cluster resemble each other, which
  # inflates the variance of the direct estimate relative to simple random
  # sampling. Simulated here by drawing contiguous blocks from the (arbitrarily
  # ordered) population frame.
  starts <- sample(seq_len(max(1, nrow(pop_i) - CLUSTER_SIZE)),
                   areas$n_clusters[i], replace = FALSE)
  idx <- unique(as.vector(sapply(starts, function(s)
    s:min(s + CLUSTER_SIZE - 1, nrow(pop_i)))))
  idx <- idx[seq_len(min(length(idx), n_take))]

  out <- pop_i[idx, ]
  out$cluster_id <- rep(seq_along(starts), each = CLUSTER_SIZE)[seq_len(nrow(out))]
  out
}

sample_parts <- lapply(seq_len(N_LGAS), sample_lga)
sample_df <- do.call(rbind, sample_parts[!vapply(sample_parts, is.null, logical(1))])
rownames(sample_df) <- NULL

# ---- 6. Population covariate means -------------------------------------------
#
# Stands in for the census or register the unit-level model requires. In real
# work this is the hard part to obtain, and the covariates must be measured
# comparably in both sources.

pop_means <- aggregate(
  cbind(hh_size, displaced, female_headed, dependency) ~ lga_id,
  population, mean
)
pop_means[, -1] <- round(pop_means[, -1], 4)

# ---- 7. Write ----------------------------------------------------------------

write.csv(areas,     file.path(OUT_DIR, "area_frame.csv"), row.names = FALSE)
write.csv(sample_df, file.path(OUT_DIR, "survey_sample.csv"), row.names = FALSE)
write.csv(pop_means, file.path(OUT_DIR, "population_covariate_means.csv"),
          row.names = FALSE)

# The full population is written so R/05 can verify against the truth. It is
# git-ignored: in a real application it does not exist, and a repository that
# shipped it would be teaching the wrong lesson.
write.csv(population, file.path(OUT_DIR, "population_TRUTH.csv"),
          row.names = FALSE)

# ---- 8. Report ---------------------------------------------------------------

cat("\nSynthetic population and sample written to", normalizePath(OUT_DIR), "\n\n")
cat("  LGAs                :", N_LGAS, "\n")
cat("  Population           :", format(nrow(population), big.mark = ","),
    "households\n")
cat("  Sample               :", format(nrow(sample_df), big.mark = ","),
    "households\n")
cat("  Sampling fraction    :",
    round(100 * nrow(sample_df) / nrow(population), 2), "%\n\n")

cat("  LGAs with no sample  :", sum(areas$n_i == 0), "of", N_LGAS, "\n")
cat("  Sample size per LGA  : min", min(areas$n_i[areas$n_i > 0]),
    " median", median(areas$n_i[areas$n_i > 0]),
    " max", max(areas$n_i), "\n")
cat("  LGAs with 0 < n < 15 :", sum(areas$n_i > 0 & areas$n_i < 15),
    "of", N_LGAS, "\n")
cat("  LGAs with 0 < n < 30 :", sum(areas$n_i > 0 & areas$n_i < 30),
    "of", N_LGAS, "\n\n")

cat("  True mean FCS        : min", round(min(areas$true_mean_fcs), 1),
    " median", round(median(areas$true_mean_fcs), 1),
    " max", round(max(areas$true_mean_fcs), 1), "\n\n")

cat("The true LGA means are retained in area_frame.csv and population_TRUTH.csv\n")
cat("so the models can be judged against them. In real work they do not exist.\n\n")
cat("None of this is real data.\n\n")
