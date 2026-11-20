# References

## Core methods

- **Fay, R.E. & Herriot, R.A. (1979).** Estimates of income for small places: an
  application of James-Stein procedures to census data. *Journal of the American
  Statistical Association* 74(366): 269-277. The area-level model implemented in
  `R/03`.

- **Battese, G.E., Harter, R.M. & Fuller, W.A. (1988).** An error-components
  model for prediction of county crop areas using survey and satellite data.
  *Journal of the American Statistical Association* 83(401): 28-36. The
  unit-level model implemented in `R/04`.

- **Rao, J.N.K. & Molina, I. (2015).** *Small Area Estimation*, 2nd edition.
  Wiley. The standard reference. Chapters 6 and 7 cover the two models and their
  MSE estimators.

- **Prasad, N.G.N. & Rao, J.N.K. (1990).** The estimation of the mean squared
  error of small-area estimators. *JASA* 85(409): 163-171. The g1 + g2 + 2*g3
  MSE decomposition used in `R/00_functions.R`.

- **Datta, G.S. & Lahiri, P. (2000).** A unified measure of uncertainty of
  estimated best linear unbiased predictors in small area estimation problems.
  *Statistica Sinica* 10: 613-627. The REML variant of the MSE.

- **Molina, I. & Rao, J.N.K. (2010).** Small area estimation of poverty
  indicators. *Canadian Journal of Statistics* 38(3): 369-385. The method to
  use for binary and non-linear indicators, which this repository does not
  implement - see the limitations section of the README.

## Validation

- **Brown, G., Chambers, R., Heady, P. & Heasman, D. (2001).** Evaluation of
  small area estimation methods: an application to unemployment estimates from
  the UK LFS. *Proceedings of Statistics Canada Symposium*. The design-based
  bias diagnostic in `R/05`.

- **Chandra, H., Salvati, N. & Chambers, R. (2018).** Small area estimation
  under a spatially non-linear model. *Computational Statistics and Data
  Analysis* 126: 19-38. On what to do when the linking model does not hold
  uniformly across space.

## Software

- **Molina, I. & Marhuenda, Y. (2015).** sae: an R package for small area
  estimation. *The R Journal* 7(1): 81-98. The production implementation. This
  repository implements the models directly so the mechanics are visible; for
  real work, use the package.

- **Bates, D., Maechler, M., Bolker, B. & Walker, S. (2015).** Fitting linear
  mixed-effects models using lme4. *Journal of Statistical Software* 67(1):
  1-48. Used for the nested error model in `R/04`.

- **Tzavidis, N. et al. (2018).** From start to finish: a framework for the
  production of small area official statistics. *Journal of the Royal
  Statistical Society A* 181(4): 927-979. A practical account of the whole
  workflow, including the reporting decisions in `R/05`.

## Survey design

- **Kish, L. (1965).** *Survey Sampling*. Wiley. The design effect used to
  correct the direct variances in `R/02`.

- **Statistics Canada (2010).** *Survey Methods and Practices*. Catalogue
  12-587-X. Source of the conventional CV publication thresholds; the specific
  values are a policy choice and are set in `R/01_config.R` rather than assumed.

## Auxiliary data in humanitarian settings

Not implemented here, but the covariates the area-level model would use in a
real application:

- VIIRS nighttime lights (NOAA/NASA), monthly composites
- WorldPop gridded population estimates
- ACLED conflict event data
- MODIS or Sentinel-derived vegetation indices
- Travel-time-to-cities accessibility surfaces (Weiss et al. 2018, *Nature* 553)
