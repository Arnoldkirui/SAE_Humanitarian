# Methodology notes

Detail that does not belong in the README but that a reviewer will ask about.

## Why the direct variance carries a design effect

Households in the same cluster resemble each other, so a clustered sample
carries less information than its raw size suggests. The design effect
quantifies it:

    deff = 1 + (b - 1) * rho

with b the average cluster size and rho the intra-cluster correlation. In the
demonstration data rho is around 0.38 at a cluster size of 6, giving a design
effect near 2.9 - the 1,854 households carry roughly the information of 640
independently drawn ones.

Omitting this is the single most common way a humanitarian survey overstates
its own precision, and it matters here for a specific reason: an understated
direct variance makes the direct estimate look reliable, which removes the
apparent need for any modelling. The comparison would then be rigged in favour
of doing nothing.

## Why REML rather than ML for the variance component

Maximum likelihood is biased downward for variance components. With 40-60
areas that bias is not negligible, and a downward-biased sigma_v2 pulls every
estimate too hard toward the regression line. The result is a set of estimates
that look impressively precise and are systematically over-shrunk - exactly
the failure a sceptical reviewer will look for.

## What a sigma_v2 of zero means

It is a real result, not a numerical failure. It says the auxiliary variables
explain the between-area variation entirely, so there is no area-specific
departure left to estimate. Every EBLUP then collapses onto the regression
line and the survey contributes nothing.

In practice it almost always means one of the auxiliaries is derived from the
outcome. `R/03` checks for this by reporting the auxiliary R-squared and
flagging values above 0.95.

## The variance of a synthetic estimate

For an area with no sample, the prediction is x'beta and its variance has two
parts:

1. the uncertainty in beta, which shrinks as the number of sampled areas grows
2. sigma_v2 itself, because that area's own departure is entirely unknown

The second part does not shrink with more data from elsewhere. In the
demonstration it is about 71% of the total variance of a synthetic estimate.

Omitting it is the standard way synthetic estimates get published looking more
precise than they are, and it is worth checking any published small-area map
for: if the unsampled areas have intervals comparable to the sampled ones,
sigma_v2 has probably been left out.

## Why a linear model is wrong for a proportion

Both models here are linear with normal errors. For a binary indicator - the
proportion of households below an FCS threshold, which is what humanitarian
reporting more often asks for - that specification fails in three ways:

1. predictions can fall outside [0,1]
2. the variance of a binary outcome depends on its mean, so homoscedastic
   errors are wrong by construction
3. the MSE estimators assume normality and do not hold

The correct approach is the empirical best prediction method of Molina & Rao
(2010), which handles non-linear indicators by simulating from the fitted
model. It is a substantially larger piece of work and this repository does not
attempt it.

Running the linear model on a proportion and reporting the result is the wrong
answer, and it is a common one.

## The sampling variance is treated as known

The Fay-Herriot model requires psi_i, the sampling variance of each direct
estimate, as a known input. It is not known - it is estimated from the same
small sample that produced the unreliable direct estimate in the first place.

In an area with six households that estimate is extremely noisy, and the model
takes it at face value when computing the shrinkage weight. A spuriously small
psi_i gives an area more weight on its own direct estimate than it deserves.

Standard practice is to smooth the variance estimates across areas before
feeding them in - fitting a model for psi_i as a function of sample size, for
example. This repository does not do that, and it is the most defensible
criticism of the implementation.

## Why the reporting rule does not use the model everywhere

The unit-level model had the lowest average RMSE in the demonstration. It would
be easy to conclude that it should be used for every area.

The reporting rule in `R/05` deliberately does not do this. Where the direct
estimate meets the CV threshold, it is published instead - because it carries
no model assumptions. When a government partner disputes a figure, "this is
what the households we interviewed reported" is a different kind of answer
from "this is what our model predicts given satellite covariates", and the
difference matters more than a small gain in average precision.

The model earns its place where the survey has nothing to offer. That is a
narrower claim than SAE is often sold on, and it is one that survives contact
with a sceptical audience.
