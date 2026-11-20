# ==============================================================================
# 01_config.R
#
# Settings for the small area estimation pipeline. Every other script sources
# this first.
#
# The thresholds below are policy choices, not statistical facts. They are
# collected here, named, so that a reviewer can disagree with a specific number
# rather than with the approach in general.
# ==============================================================================

# ------------------------------------------------------------------------------
# EDIT THIS BLOCK
# ------------------------------------------------------------------------------

DATA_DIR   <- "data"
OUTPUT_DIR <- "output"

# The indicator being estimated, and the area it is estimated for
INDICATOR  <- "Mean Food Consumption Score"
AREA_LABEL <- "LGA"

# --- Reporting thresholds -----------------------------------------------------
#
# Coefficient of variation above which an estimate is flagged or withheld.
# Statistical offices commonly publish below 25%, flag between 25 and 35, and
# suppress above 35. There is no universal rule; what matters is that the
# threshold is set in advance and stated in the methodology, rather than chosen
# after seeing which areas it excludes.

CV_PUBLISH_THRESHOLD <- 25    # at or below: publish
CV_FLAG_THRESHOLD    <- 35    # above this: suppress

# --- Auxiliary variables for the area-level model -----------------------------
#
# These must be available for EVERY area, including the ones with no survey
# sample at all. That requirement rules out anything derived from the survey
# itself, which is the most common way an area-level model is accidentally
# invalidated.

FH_COVARIATES <- c("log_nightlights", "log_pop_density",
                   "market_distance", "sqrt_conflict", "ndvi_anomaly")

# --- Covariates for the unit-level model --------------------------------------
#
# Each of these needs a POPULATION mean per area, from a census or register.
# The survey cannot supply it - if it could, the area would not be small.

BHF_COVARIATES <- c("hh_size", "displaced", "female_headed", "dependency")

# --- Validation ---------------------------------------------------------------

# Set FALSE to run the pipeline as it would run on real data, where the true
# area values do not exist and R/05's comparison against them is impossible.
TRUTH_AVAILABLE <- TRUE

# ------------------------------------------------------------------------------
# NOTHING BELOW NEEDS EDITING
# ------------------------------------------------------------------------------

stopifnot(CV_FLAG_THRESHOLD > CV_PUBLISH_THRESHOLD)

for (d in c(DATA_DIR, OUTPUT_DIR,
            file.path(OUTPUT_DIR, "tables"),
            file.path(OUTPUT_DIR, "figures"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

message("Indicator      : ", INDICATOR)
message("Estimated for  : ", AREA_LABEL)
message("CV thresholds  : publish <= ", CV_PUBLISH_THRESHOLD,
        "%, suppress > ", CV_FLAG_THRESHOLD, "%")
if (TRUTH_AVAILABLE) {
  message("\nSYNTHETIC DEMONSTRATION DATA - the true area values are known,")
  message("which is why the models can be validated. On real data they are not.")
}
