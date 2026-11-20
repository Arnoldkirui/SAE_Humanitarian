# Small Area Estimation for Humanitarian Assessments

Household surveys are powered for the level they report at. A national MSNA is
designed to describe states, and then everyone asks for LGA figures.

The honest answer has usually been that the cell is too small to report. Small
area estimation is the method that turns some of those cells into defensible
numbers, and tells you which ones it cannot.

This repository implements both standard approaches, applies them to a
simulated assessment where the true values are known, and then checks whether
they actually worked.

**Status: methodology and code. The data is synthetic by design.**

## The problem, in the demonstration data

A survey of 1,854 households across 60 LGAs, allocated the way a real
state-level assessment allocates:

| | LGAs | Households |
|---|---|---|
| Direct estimate publishable (CV ≤ 25%) | 46 | 92,237 |
| Direct estimate too imprecise | 5 | 6,991 |
| **No sample at all** | **9** | **16,594** |

Fourteen of sixty LGAs have no usable figure, covering about a fifth of the
population. The assessment happened; the numbers cannot carry weight.

## What this does

| Script | What it produces |
|---|---|
| `setup_data/00_generate_population.R` | A synthetic finite population with known true LGA means, plus a realistic two-stage sample drawn from it |
| `R/02_direct_estimates.R` | Direct estimates with design-effect-corrected variances, and the CV screen that establishes the problem |
| `R/03_fay_herriot.R` | Area-level model: borrows strength from satellite and administrative auxiliaries |
| `R/04_unit_level_bhf.R` | Unit-level model: nested error model on household records |
| `R/05_compare_and_validate.R` | Both models judged against the truth, plus the diagnostics that work without it |

## The finding

Running both and checking against the truth produces a result that an SAE
demonstration does not usually report:

| Subset | Direct | Fay–Herriot | Unit-level |
|---|---|---|---|
| All sampled areas | **2.37** | 2.40 | 2.32 |
| Sampled, n < 15 | 3.68 | 3.56 | **3.52** |
| Sampled, n ≥ 30 | **1.55** | 1.73 | 1.59 |
| Unsampled areas | — | **4.59** | 8.06 |

*RMSE against the true LGA mean, FCS points.*

Across sampled areas as a whole, **the area-level model is not more accurate
than the direct estimate.** Shrinkage helps the thin areas and costs something
in the well-sampled ones, and the average comes out flat.

That is not a failure of the implementation, and it matters because it changes
what can honestly be claimed. The case for SAE is not "it beats the survey". It
is narrower: it rescues areas the survey cannot support, and it should leave the
rest alone. The reporting rule in `R/05` is built on that narrower claim.

The unsampled row is the other finding. The two models are far apart, and the
reason is covariates rather than model class. The area-level model uses
auxiliaries that vary between areas — nightlights, conflict incidents, market
access — so it can say something specific about a place it has never visited.
The unit-level model uses household characteristics, which predict a
household's outcome well but carry little information about why one LGA differs
from another; its synthetic prediction collapses toward the national average and
picks up a 4.2-point bias.

For areas with no sample, what matters is whether your covariates explain
**between-area** variation. A more sophisticated model on covariates that do not
will not help.

## Where the model makes things worse

22 of 51 sampled areas have a Fay–Herriot estimate further from the truth than
the direct one. `R/05` lists the worst and does not explain them away.

Shrinkage trades a worse estimate in genuinely unusual areas for better
estimates on average — and an area that truly departs from the pattern is
exactly the one it damages. If the outlying LGAs are the ones a programme cares
about, and in humanitarian work the worst-off usually are, that trade may not be
acceptable. The right answer is then to sample those areas properly rather than
to model them.

## Validating without the truth

Everything above needs the true values, which exist only because the population
is simulated. On real data the model's claim about its own accuracy is internal
— an MSE derived from the model's own assumptions, which is not evidence that
the assumptions hold.

`R/05` therefore also runs the diagnostics that transfer:

- **Design-based bias check.** In areas where the direct estimate *is* reliable,
  the model should agree with it. Systematic disagreement where the survey can
  be trusted is a reason to stop.
- **Interval coverage.** Standardised differences between the two estimates,
  against what their stated variances imply.

## Getting started

```r
install.packages(c("dplyr", "readr", "tidyr", "ggplot2", "lme4"))

source("setup_data/00_generate_population.R")

source("R/02_direct_estimates.R")
source("R/03_fay_herriot.R")
source("R/04_unit_level_bhf.R")
source("R/05_compare_and_validate.R")
```

Open `SAE_Humanitarian.Rproj` first so relative paths resolve. All scripts have
been run end to end under R 4.3.3.

## Choosing between the two models

**Area-level (Fay–Herriot)** needs only the direct estimates, their sampling
variances, and auxiliary data per area. It is usually the right choice in
humanitarian settings, because satellite and administrative covariates are
available where census microdata is not.

**Unit-level (Battese–Harter–Fuller)** is more precise when it applies, but
needs the **population mean of every covariate in every area** — from a census
or register, measured comparably to your survey. A roster counted one way in
your questionnaire and another in the census is not the same variable, and using
it as though it were introduces a bias no diagnostic here will catch.

That requirement, not the mathematics, is why area-level models dominate this
sector.

## Limitations

- **The models are linear and assume normal errors.** That suits the continuous
  FCS used here. It would **not** suit a binary indicator such as the proportion
  below a threshold, which is what humanitarian reporting more often wants. That
  needs a generalised linear mixed model and the empirical best prediction
  literature, which this repository does not cover.
- **Sampling variances are treated as known**, as the Fay–Herriot model
  requires. In a small area they are estimated from very few observations and
  are themselves uncertain; the model does not propagate that.
- **The unit-level MSE carries leading terms only.** A bootstrap is the honest
  route for a real application.
- **Not validated on real data.** The true values here are known by
  construction. Nothing in this repository establishes how either model behaves
  on an actual assessment.

## Licence

Code MIT. Documentation CC BY 4.0.

## Author

[Arnold Kirui](https://github.com/Arnoldkirui), a data specialist with
experience across East and West Africa — Uganda, South Sudan, Kenya and
Nigeria — working in the humanitarian sector.

Corrections welcome, particularly from anyone who has run SAE against a real
assessment and can say how the design-based checks behaved.
