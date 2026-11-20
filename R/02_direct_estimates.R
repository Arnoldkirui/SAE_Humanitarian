# ==============================================================================
# 02_direct_estimates.R
#
# Computes the direct survey estimate for every LGA, with a variance that
# accounts for clustering. These are the estimates a standard analysis would
# produce and the baseline everything else is judged against.
#
# This script's job is to establish the problem, not solve it. The output to
# look at is the CV column: in the thin LGAs it is large enough that the
# estimate cannot responsibly be published, which is the situation that makes
# small area estimation worth doing.
#
# Run after setup_data/00_generate_population.R.
# ==============================================================================

source("R/01_config.R")
source("R/00_functions.R")

library(dplyr)
library(readr)
library(ggplot2)


# ==============================================================================
# 1. Load
# ==============================================================================

sample_df <- read_csv(file.path(DATA_DIR, "survey_sample.csv"),
                      show_col_types = FALSE)
areas     <- read_csv(file.path(DATA_DIR, "area_frame.csv"),
                      show_col_types = FALSE)

n_sampled   <- n_distinct(sample_df$lga_id)
n_unsampled <- nrow(areas) - n_sampled

message("\nSample: ", nrow(sample_df), " households in ",
        n_sampled, " of ", nrow(areas), " LGAs")

if (n_unsampled > 0) {
  no_sample <- areas |> filter(!lga_id %in% unique(sample_df$lga_id))
  message("\n", n_unsampled, " LGA(s) have NO sample at all, covering ",
          format(sum(no_sample$N_i), big.mark = ","), " households.")
  message("They have no direct estimate - not an imprecise one, none. ",
          "Nothing in this script can say anything about them.")
}


# ==============================================================================
# 2. Intra-cluster correlation
# ==============================================================================
#
# Households in the same cluster resemble each other, so each interview carries
# less independent information than a simple random sample would. Ignoring this
# understates the variance of the direct estimate - and then the direct
# estimate looks more reliable than it is, which removes the apparent need for
# any of what follows.
#
# ICC is estimated from a one-way random effects decomposition across all
# clusters, pooled. Estimating it per LGA would be preferable and is not
# possible here: most LGAs have too few clusters.

icc_fit <- aov(fcs ~ factor(paste(lga_id, cluster_id)), data = sample_df)
ms <- summary(icc_fit)[[1]][["Mean Sq"]]

k <- mean(table(paste(sample_df$lga_id, sample_df$cluster_id)))
icc <- max(0, (ms[1] - ms[2]) / (ms[1] + (k - 1) * ms[2]))

deff <- design_effect(k, icc)

message("\nIntra-cluster correlation: ", round(icc, 4))
message("Average cluster size     : ", round(k, 1))
message("Design effect            : ", round(deff, 3))
message("\nThat design effect means the sample carries about the information ",
        "of ", round(nrow(sample_df) / deff), " independently drawn ",
        "households, not ", nrow(sample_df), ".")


# ==============================================================================
# 3. Direct estimates
# ==============================================================================
#
# The sample is self-weighting within LGA here, so the direct estimator is the
# sample mean. A real MSNA carries design weights, and the direct estimate
# would be the weighted mean with a variance from the design - the structure
# below is unchanged, only the two formulas differ.

direct <- sample_df |>
  group_by(lga_id) |>
  summarise(
    n_i        = n(),
    n_clusters = n_distinct(cluster_id),
    direct_est = mean(fcs),
    s2         = var(fcs),
    .groups = "drop"
  ) |>
  mutate(
    # Variance of the mean, inflated by the design effect
    var_srs    = s2 / n_i,
    direct_var = var_srs * deff,
    direct_se  = sqrt(direct_var),
    cv         = cv_pct(direct_est, direct_var)
  )

# An LGA with a single cluster has no between-cluster variation to estimate
# from, so its variance is understated by construction. Flag it rather than
# letting it quietly look precise.
single_cluster <- direct |> filter(n_clusters < 2)
if (nrow(single_cluster) > 0) {
  message("\n", nrow(single_cluster), " LGA(s) have fewer than 2 clusters. ",
          "Their variance is not properly estimable from the design and is ",
          "understated here.")
}

direct <- direct |>
  left_join(areas |> select(lga_id, lga_name, state, N_i, true_mean_fcs),
            by = "lga_id")


# ==============================================================================
# 4. The problem, stated in numbers
# ==============================================================================

reportable <- direct |>
  mutate(
    status = case_when(
      cv <= CV_PUBLISH_THRESHOLD  ~ "publishable",
      cv <= CV_FLAG_THRESHOLD     ~ "flag with caution",
      TRUE                        ~ "suppress"
    )
  )

message("\n=============================================================")
message("DIRECT ESTIMATES AGAINST THE REPORTING THRESHOLDS")
message("=============================================================\n")
print(as.data.frame(count(reportable, status)))

message("\nCV distribution: min ", round(min(direct$cv), 1),
        "%, median ", round(median(direct$cv), 1),
        "%, max ", round(max(direct$cv), 1), "%")

n_unusable <- sum(reportable$status != "publishable")
hh_unusable <- sum(reportable$N_i[reportable$status != "publishable"])
if (n_unsampled > 0) hh_unusable <- hh_unusable + sum(no_sample$N_i)

message("\n", n_unusable, " of ", nrow(direct), " sampled LGAs cannot be ",
        "published from the direct estimate at a ", CV_PUBLISH_THRESHOLD,
        "% CV threshold.")
message("Adding the ", n_unsampled, " unsampled LGAs, ",
        n_unusable + n_unsampled, " of ", nrow(areas),
        " LGAs have no usable figure - ",
        format(hh_unusable, big.mark = ","), " households.")

message("\nThat is the operational problem, and it has two halves. In the ",
        "thin LGAs the data exists and the estimate is too imprecise to act ",
        "on. In the unsampled ones there is nothing at all.")
message("The options are to publish figures that cannot carry the weight, ",
        "to say nothing about a quarter of the population, or to model.")


# ==============================================================================
# 5. How wrong are the direct estimates?
# ==============================================================================
#
# Only answerable because this is simulated data. On a real survey this section
# cannot exist, which is precisely why the demonstration is built this way.

direct <- direct |>
  mutate(
    direct_error = direct_est - true_mean_fcs,
    abs_error    = abs(direct_error)
  )

error_by_size <- direct |>
  mutate(size_band = cut(n_i, breaks = c(0, 12, 24, 42, Inf),
                         labels = c("6-12", "13-24", "25-42", "43+"))) |>
  group_by(size_band) |>
  summarise(
    n_lgas        = n(),
    mean_abs_error = round(mean(abs_error), 2),
    max_abs_error  = round(max(abs_error), 2),
    mean_cv        = round(mean(cv), 1),
    .groups = "drop"
  )

message("\nAbsolute error of the direct estimate, by sample size:")
print(as.data.frame(error_by_size))

message("\nThe error falls steadily with sample size, as it must. The point of ",
        "the table is the top row: in the thinnest LGAs the direct estimate ",
        "misses the truth by several FCS points, which is enough to move an ",
        "area across a programming threshold.")


# ==============================================================================
# 6. Outputs
# ==============================================================================

write_csv(direct, file.path(OUTPUT_DIR, "tables", "direct_estimates.csv"))
write_csv(error_by_size,
          file.path(OUTPUT_DIR, "tables", "direct_error_by_size.csv"))

saveRDS(list(icc = icc, deff = deff, avg_cluster = k),
        file.path(OUTPUT_DIR, "design_parameters.rds"))

p <- ggplot(direct, aes(x = n_i, y = cv)) +
  geom_point(colour = "#2E6E8E", size = 2, alpha = 0.85) +
  geom_hline(yintercept = CV_PUBLISH_THRESHOLD, linetype = "dashed",
             colour = "#B2182B") +
  geom_hline(yintercept = CV_FLAG_THRESHOLD, linetype = "dotted",
             colour = "grey40") +
  annotate("text", x = max(direct$n_i) * 0.78,
           y = CV_PUBLISH_THRESHOLD + 1.6,
           label = paste0(CV_PUBLISH_THRESHOLD, "% publication threshold"),
           size = 3, colour = "#B2182B") +
  labs(
    title = "Most LGAs cannot be reported from the direct estimate",
    subtitle = paste0("Coefficient of variation against LGA sample size. ",
                      n_unusable, " of ", nrow(direct),
                      " LGAs sit above the threshold."),
    x = "Households interviewed in the LGA",
    y = "CV of the direct estimate (%)",
    caption = paste0("Design effect ", round(deff, 2),
                     " applied. Synthetic demonstration data.")
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title.position = "plot",
        plot.caption = element_text(hjust = 0, colour = "grey45"))

ggsave(file.path(OUTPUT_DIR, "figures", "direct_cv_by_sample_size.png"),
       p, width = 20, height = 12, units = "cm", dpi = 300)

message("\nWritten to ", OUTPUT_DIR, "/tables and ", OUTPUT_DIR, "/figures")
