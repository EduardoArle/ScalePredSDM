# ScalePredSDM

ScalePredSDM is an R package providing a reproducible workflow for generating,
evaluating and applying scale-of-effect (SoE) predictors in fine-scale species
distribution models (SDMs).

The package is designed for ecological applications where predictor scale
selection is critical, including habitat suitability modelling, seascape
ecology and analyses involving positional uncertainty in occurrence data.

---

## Scientific background

Ecological responses to environmental predictors are often scale-dependent.
Conventional SDM workflows typically rely on a single, user-defined spatial
resolution, which can introduce bias and reduce transferability.  
ScalePredSDM supports a systematic, data-driven approach to estimating
species-specific scales of effect while maintaining reproducibility and
compatibility with modern SDM workflows.

The package was developed with a focus on applied SDMs, where methodological transparency and
consistency are essential.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("olejohs/ScalePredSDM")
