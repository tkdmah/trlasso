# Least Absolute Shrinkage and Selection Operator with Transfer Learning

[(Takada & Fujisawa, 2020, NeurIPS)](https://papers.nips.cc/paper/2020/hash/a4a83056b58ff983d12c72bb17996243-Abstract.html) proposed a method for transferring knowledge from a source domain to a target domain via $\ell_1$ regularization in high dimension. [(Takada & Fujisawa 2023 arXiv)](https://arxiv.org/abs/2308.15838) extended the method to the weighted regularization to improve its theoretical properties. This is a simple implementation of the above methods.

## Installation
To install the package:
```
# install.packages("remotes")
remotes::install_github("tkdmah/trlasso", build_vignettes = TRUE, dependencies = TRUE)
```

## Usage
To learn how to use:
```
vignette("example_linear", package = "trlasso")
vignette("example_linear_weight", package = "trlasso")
```

| File | Description |
| --- | --- |
|vignettes/example_linear.Rmd | Sample script for Transfer Lasso|
|vignettes/example_linear_weight.Rmd | Sample script for Weighted Transfer Lasso|

## Experiment Scripts
To replicate the experiments in [(Takada & Fujisawa, 2023, arXiv)](https://arxiv.org/abs/2308.15838):
| File | Description |
| --- | --- |
| inst/experiments/empirical_convergence.Rmd | Experiment script for evaluating estimation errors with respect to sample size |
| inst/experiments/empirical_phase_diagram.Rmd | Experiment script for drawing empirical phase diagrams |
| inst/experiments/method_comparison_large.Rmd | Experiment script for comparing several methods using a large amount of source data (Note: time-consuming) |
| inst/experiments/method_comparison_medium.Rmd | Experiment script for comparing several methods using a medium amount of source data (Note: time-consuming) |

## References
1. [Takada, M., & Fujisawa, H. (2020). Transfer Learning via $\ell_1$ Regularization. Advances in Neural Information Processing Systems, 33, 14266-14277.](https://papers.nips.cc/paper/2020/hash/a4a83056b58ff983d12c72bb17996243-Abstract.html)
2. [Takada, M., & Fujisawa, H. (2023). Adaptive Lasso, Transfer Lasso, and Beyond: An Asymptotic Perspective. arXiv preprint arXiv:2308.15838.](https://arxiv.org/abs/2308.15838)

## Developers
- Masaaki Takada (Toshiba Corporation)
- Gen Li (Toshiba Corporation)
- Sunao Yotsutsuji (Toshiba Corporation)
