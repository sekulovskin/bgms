# Extract Coefficients from a bgms Object

Returns the posterior mean main effects, pairwise effects, and edge
inclusion indicators from a `bgms` model fit.

## Usage

``` r
# S3 method for class 'bgms'
coef(object, ...)
```

## Arguments

- object:

  An object of class `bgms`.

- ...:

  Ignored.

## Value

A list with the following components:

- main:

  Posterior mean of the main-effect parameters. `NULL` for GGM models
  (no main effects). For OMRF models this is a numeric matrix (p x
  max_categories) of category thresholds. For mixed MRF models this is a
  list with `$discrete` (p x max_categories matrix) and `$continuous` (q
  x 1 matrix of means).

- pairwise:

  Posterior mean of the pairwise interaction matrix. For GGM and mixed
  MRF models the precision matrix diagonal is included on the matrix
  diagonal.

- indicator:

  Posterior mean of the edge inclusion indicators (if available).

## See also

[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
[`print.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

Other posterior-methods:
[`coef.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgmCompare.md),
[`print.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md),
[`print.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md),
[`summary.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

## Examples

``` r
# \donttest{
fit = bgm(x = Wenchuan[, 1:3])
#> 2 rows with missing values excluded (n = 360 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/2000 (7.5%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 170/2000 (8.5%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 167/2000 (8.3%)
#> Chain 4 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 224/2000 (11.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 711/8000 (8.9%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 519/2000 (25.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 513/2000 (25.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 579/2000 (28.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2061/8000 (25.8%)
#> Elapsed: 1s | ETA: 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 882/2000 (44.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 880/2000 (44.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 934/2000 (46.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3496/8000 (43.7%)
#> Elapsed: 1s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1197/2000 (59.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1192/2000 (59.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1256/2000 (62.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4745/8000 (59.3%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1450/2000 (72.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1553/2000 (77.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1535/2000 (76.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1609/2000 (80.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6147/8000 (76.8%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1899/2000 (95.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 1856/2000 (92.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 1952/2000 (97.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7507/8000 (93.8%)
#> Elapsed: 3s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 3s | ETA: 0s
coef(fit)
#> $main
#>               cat (1)   cat (2)   cat (3)    cat (4)
#> intrusion  0.80908860 -1.172276 -3.530161  -7.430128
#> dreams    -0.39966517 -3.299039 -6.365205 -10.468118
#> flash      0.07203052 -2.127073 -4.573521  -8.406048
#> 
#> $pairwise
#>           intrusion    dreams     flash
#> intrusion 0.0000000 0.7010837 0.4236871
#> dreams    0.7010837 0.0000000 0.5617632
#> flash     0.4236871 0.5617632 0.0000000
#> 
#> $indicator
#>           intrusion dreams flash
#> intrusion         0      1     1
#> dreams            1      0     1
#> flash             1      1     0
#> 
# }
```
