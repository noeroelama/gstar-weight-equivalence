# Extra analyses (generated 2026-10-04 04:48)

## A/H. Five candidates incl. uniform W; dynamic-part margin

| W                    |   phi11_mean |   direct |   indirect |   rho_A |
|:---------------------|-------------:|---------:|-----------:|--------:|
| geographic KNN       |        0.161 |    1.155 |      0.214 |   0.456 |
| inverse-distance     |        0.243 |    1.135 |      0.397 |   0.401 |
| training-correlation |        0.265 |    1.134 |      0.474 |   0.401 |
| graphical-lasso      |        0.196 |    1.123 |      0.257 |   0.487 |
| uniform              |        0.243 |    1.166 |      0.444 |   0.4   |

|   theta_hat |   se_mbb |   L_pw |   ci90_lo |   ci90_hi |   z_se |   tau_ref |   c_min | equivalent   |   c_min_a025 | equivalent_a025   | equivalent_hac   |   dyn_tau_ref |   dyn_c_min |   dyn_equivalent |   dyn_c_min_a025 |   dyn_equivalent_a025 |   dyn_equivalent_hac | functional   | W_A                  | W_B                  |
|------------:|---------:|-------:|----------:|----------:|-------:|----------:|--------:|:-------------|-------------:|:------------------|:-----------------|--------------:|------------:|-----------------:|-----------------:|----------------------:|---------------------:|:-------------|:---------------------|:---------------------|
|       0.02  |    0.014 |  2.748 |    -0.003 |     0.043 |  0.023 |     1.155 |   0.038 | True         |        0.041 | True              | True             |         0.155 |       0.271 |                0 |            0.299 |                     0 |                    0 | direct       | geographic KNN       | inverse-distance     |
|       0.02  |    0.016 |  2.748 |    -0.005 |     0.046 |  0.026 |     1.155 |   0.04  | True         |        0.044 | True              | True             |         0.155 |       0.276 |                0 |            0.302 |                     0 |                    0 | direct       | geographic KNN       | training-correlation |
|       0.032 |    0.011 |  2.748 |     0.014 |     0.05  |  0.018 |     1.155 |   0.043 | True         |        0.046 | True              | True             |         0.155 |       0.302 |                0 |            0.32  |                     0 |                    0 | direct       | geographic KNN       | graphical-lasso      |
|      -0.011 |    0.016 |  2.748 |    -0.037 |     0.015 |  0.026 |     1.155 |   0.033 | True         |        0.038 | True              | True             |         0.155 |       0.278 |                0 |            0.326 |                     0 |                    0 | direct       | geographic KNN       | uniform              |
|       0.001 |    0.006 |  2.748 |    -0.009 |     0.01  |  0.01  |     1.155 |   0.009 | True         |        0.01  | True              | True             |         0.155 |       0.068 |                1 |            0.084 |                     1 |                    1 | direct       | inverse-distance     | training-correlation |
|       0.012 |    0.013 |  2.964 |    -0.009 |     0.033 |  0.021 |     1.155 |   0.029 | True         |        0.032 | True              | True             |         0.155 |       0.22  |                0 |            0.25  |                     0 |                    0 | direct       | inverse-distance     | graphical-lasso      |
|      -0.031 |    0.01  |  2.748 |    -0.048 |    -0.015 |  0.016 |     1.155 |   0.042 | True         |        0.044 | True              | True             |         0.155 |       0.336 |                0 |            0.366 |                     0 |                    0 | direct       | inverse-distance     | uniform              |
|       0.011 |    0.013 |  2.748 |    -0.01  |     0.033 |  0.021 |     1.155 |   0.028 | True         |        0.032 | True              | True             |         0.155 |       0.222 |                0 |            0.256 |                     0 |                    0 | direct       | training-correlation | graphical-lasso      |
|      -0.032 |    0.007 |  2.748 |    -0.043 |    -0.021 |  0.011 |     1.155 |   0.037 | True         |        0.038 | True              | True             |         0.155 |       0.288 |                0 |            0.308 |                     0 |                    0 | direct       | training-correlation | uniform              |
|      -0.043 |    0.015 |  2.748 |    -0.068 |    -0.018 |  0.025 |     1.155 |   0.059 | True         |        0.064 | True              | True             |         0.155 |       0.475 |                0 |            0.52  |                     0 |                    0 | direct       | graphical-lasso      | uniform              |
|      -0.183 |    0.069 |  3.535 |    -0.297 |    -0.069 |  0.114 |     0.214 |   1.233 | False        |        1.314 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | inverse-distance     |
|      -0.26  |    0.088 |  3.148 |    -0.404 |    -0.116 |  0.144 |     0.214 |   1.707 | False        |        1.821 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | training-correlation |
|      -0.043 |    0.017 |  1.388 |    -0.072 |    -0.014 |  0.029 |     0.214 |   0.381 | False        |        0.432 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | graphical-lasso      |
|      -0.23  |    0.088 |  3.391 |    -0.375 |    -0.085 |  0.145 |     0.214 |   1.637 | False        |        1.769 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | uniform              |
|      -0.076 |    0.023 |  1.388 |    -0.114 |    -0.039 |  0.038 |     0.214 |   0.57  | False        |        0.63  | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | inverse-distance     | training-correlation |
|       0.14  |    0.061 |  3.75  |     0.039 |     0.241 |  0.101 |     0.214 |   0.992 | False        |        1.063 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | inverse-distance     | graphical-lasso      |
|      -0.047 |    0.028 |  2.818 |    -0.093 |    -0.001 |  0.046 |     0.214 |   0.473 | False        |        0.538 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | inverse-distance     | uniform              |
|       0.216 |    0.08  |  3.377 |     0.086 |     0.347 |  0.131 |     0.214 |   1.438 | False        |        1.531 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | training-correlation | graphical-lasso      |
|       0.03  |    0.016 |  1.388 |     0.003 |     0.056 |  0.026 |     0.214 |   0.263 | False        |        0.294 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | training-correlation | uniform              |
|      -0.187 |    0.081 |  3.613 |    -0.32  |    -0.054 |  0.133 |     0.214 |   1.37  | False        |        1.48  | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | graphical-lasso      | uniform              |

## B. Common-factor augmentation (lagged cross-sectional mean)

| W                    |   phi11_mean |   phi_mean_coef |   direct |   indirect |   rho_A |   indirect_W_only |
|:---------------------|-------------:|----------------:|---------:|-----------:|--------:|------------------:|
| geographic KNN       |        0.052 |           0.211 |    1.154 |      0.445 |   0.478 |             0.031 |
| inverse-distance     |       -0.034 |           0.272 |    1.195 |      0.419 |   0.518 |            -0.125 |
| training-correlation |        0.526 |          -0.258 |    1.128 |      0.471 |   0.417 |             1.332 |
| graphical-lasso      |        0.077 |           0.196 |    1.145 |      0.459 |   0.497 |             0.058 |

|   theta_hat |   se_mbb |   L_pw |   ci90_lo |   ci90_hi |   z_se |   tau_ref |   c_min | equivalent   |   c_min_a025 | equivalent_a025   | equivalent_hac   |   dyn_tau_ref |   dyn_c_min |   dyn_equivalent |   dyn_c_min_a025 |   dyn_equivalent_a025 |   dyn_equivalent_hac | functional   | W_A                  | W_B                  |
|------------:|---------:|-------:|----------:|----------:|-------:|----------:|--------:|:-------------|-------------:|:------------------|:-----------------|--------------:|------------:|-----------------:|-----------------:|----------------------:|---------------------:|:-------------|:---------------------|:---------------------|
|      -0.041 |    0.013 |  2.658 |    -0.062 |    -0.019 |  0.021 |     1.154 |   0.054 | True         |        0.057 | True              | True             |         0.154 |       0.4   |                0 |            0.43  |                     0 |                    0 | direct       | geographic KNN       | inverse-distance     |
|       0.026 |    0.015 |  2.658 |     0.001 |     0.05  |  0.025 |     1.154 |   0.044 | True         |        0.047 | True              | True             |         0.154 |       0.313 |                0 |            0.342 |                     0 |                    0 | direct       | geographic KNN       | training-correlation |
|       0.009 |    0.011 |  2.658 |    -0.009 |     0.026 |  0.017 |     1.154 |   0.022 | True         |        0.026 | True              | True             |         0.154 |       0.161 |                0 |            0.18  |                     0 |                    0 | direct       | geographic KNN       | graphical-lasso      |
|       0.066 |    0.02  |  2.658 |     0.034 |     0.099 |  0.033 |     1.154 |   0.084 | True         |        0.09  | True              | True             |         0.154 |       0.602 |                0 |            0.637 |                     0 |                    0 | direct       | inverse-distance     | training-correlation |
|       0.049 |    0.018 |  2.658 |     0.02  |     0.079 |  0.029 |     1.154 |   0.067 | True         |        0.072 | True              | True             |         0.154 |       0.482 |                0 |            0.514 |                     0 |                    0 | direct       | inverse-distance     | graphical-lasso      |
|      -0.017 |    0.013 |  2.658 |    -0.038 |     0.004 |  0.021 |     1.154 |   0.034 | True         |        0.037 | True              | True             |         0.154 |       0.252 |                0 |            0.282 |                     0 |                    0 | direct       | training-correlation | graphical-lasso      |
|       0.026 |    0.013 |  2.432 |     0.005 |     0.048 |  0.022 |     0.445 |   0.142 | False        |        0.17  | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | inverse-distance     |
|      -0.026 |    0.016 |  2.432 |    -0.053 |     0.001 |  0.027 |     0.445 |   0.136 | False        |        0.158 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | training-correlation |
|      -0.014 |    0.012 |  2.432 |    -0.033 |     0.005 |  0.019 |     0.445 |   0.083 | True         |        0.097 | True              | True             |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | geographic KNN       | graphical-lasso      |
|      -0.052 |    0.019 |  2.432 |    -0.083 |    -0.022 |  0.031 |     0.445 |   0.237 | False        |        0.278 | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | inverse-distance     | training-correlation |
|      -0.04  |    0.017 |  2.432 |    -0.068 |    -0.012 |  0.028 |     0.445 |   0.188 | False        |        0.22  | False             | False            |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | inverse-distance     | graphical-lasso      |
|       0.012 |    0.015 |  2.703 |    -0.013 |     0.037 |  0.025 |     0.445 |   0.098 | True         |        0.118 | False             | True             |       nan     |     nan     |              nan |          nan     |                   nan |                  nan | indirect     | training-correlation | graphical-lasso      |

## C. Mean correlation between spatial-lag series

|                      |   geographic KNN |   inverse-distance |   training-correlation |   graphical-lasso |   uniform |
|:---------------------|-----------------:|-------------------:|-----------------------:|------------------:|----------:|
| geographic KNN       |            1     |              0.907 |                  0.836 |             0.954 |     0.743 |
| inverse-distance     |            0.907 |              1     |                  0.979 |             0.885 |     0.941 |
| training-correlation |            0.836 |              0.979 |                  1     |             0.832 |     0.981 |
| graphical-lasso      |            0.954 |              0.885 |                  0.832 |             1     |     0.735 |
| uniform              |            0.743 |              0.941 |                  0.981 |             0.735 |     1     |

## D. J-test

| maintained           | rival                |   share_reject_5pct |   mean_rho |   rho_pooled |   se_pooled_mbb |   t_pooled | reject_pooled   |
|:---------------------|:---------------------|--------------------:|-----------:|-------------:|----------------:|-----------:|:----------------|
| geographic KNN       | inverse-distance     |               0.3   |      1.25  |        0.859 |           0.239 |      3.594 | True            |
| geographic KNN       | training-correlation |               0.5   |      0.869 |        0.91  |           0.189 |      4.807 | True            |
| geographic KNN       | graphical-lasso      |               0.375 |      0.219 |        0.985 |           0.172 |      5.717 | True            |
| inverse-distance     | geographic KNN       |               0.125 |      0.049 |        0.567 |           0.154 |      3.693 | True            |
| inverse-distance     | training-correlation |               0.425 |      1.176 |        1.899 |           0.405 |      4.689 | True            |
| inverse-distance     | graphical-lasso      |               0.275 |      0.425 |        0.699 |           0.136 |      5.138 | True            |
| training-correlation | geographic KNN       |               0.15  |      0.293 |        0.581 |           0.118 |      4.929 | True            |
| training-correlation | inverse-distance     |               0.35  |     -1.945 |       -0.62  |           0.381 |     -1.625 | False           |
| training-correlation | graphical-lasso      |               0.2   |      0.164 |        0.592 |           0.111 |      5.357 | True            |
| graphical-lasso      | geographic KNN       |               0.25  |     -1.523 |        0.672 |           0.162 |      4.141 | True            |
| graphical-lasso      | inverse-distance     |               0.3   |      1.223 |        0.821 |           0.214 |      3.845 | True            |
| graphical-lasso      | training-correlation |               0.475 |      0.983 |        0.899 |           0.191 |      4.703 | True            |

## E. Benchmarks

| model                                          |   rmse_mm |   ci90_lo |   ci90_hi |
|:-----------------------------------------------|----------:|----------:|----------:|
| climatology only                               |    80.846 |   nan     |   nan     |
| AR(1) + climate, no W                          |    74.321 |   nan     |   nan     |
| GSTARX geographic KNN                          |    73.996 |   nan     |   nan     |
| GSTARX inverse-distance                        |    73.687 |   nan     |   nan     |
| GSTARX training-correlation                    |    73.419 |   nan     |   nan     |
| GSTARX graphical-lasso                         |    73.996 |   nan     |   nan     |
| GSTARX uniform                                 |    73.502 |   nan     |   nan     |
| GSTARX geographic KNN minus AR(1): theta       |    -0.325 |    -0.53  |    -0.121 |
| GSTARX inverse-distance minus AR(1): theta     |    -0.635 |    -0.978 |    -0.291 |
| GSTARX training-correlation minus AR(1): theta |    -0.903 |    -1.347 |    -0.458 |
| GSTARX graphical-lasso minus AR(1): theta      |    -0.325 |    -0.582 |    -0.068 |
| GSTARX uniform minus AR(1): theta              |    -0.82  |    -1.248 |    -0.391 |

## F. Sensitivity

| setting                        |   direct |   indirect |   rho_A |   theta_vs_corr_indirect |   c_min |   overlap_with_a_geo |   rmse_knn |   rmse_corr |   M |   theta_rmse |   min_delta |   equivalent_1mm |
|:-------------------------------|---------:|-----------:|--------:|-------------------------:|--------:|---------------------:|-----------:|------------:|----:|-------------:|------------:|-----------------:|
| KNN k=3 (centroid-based)       |    1.177 |      0.161 |   0.47  |                   -0.313 |   2.845 |                  0.6 |    nan     |     nan     | nan |      nan     |     nan     |              nan |
| KNN k=5 (centroid-based)       |    1.155 |      0.214 |   0.456 |                   -0.26  |   1.708 |                  1   |    nan     |     nan     | nan |      nan     |     nan     |              nan |
| KNN k=8 (centroid-based)       |    1.157 |      0.235 |   0.457 |                   -0.239 |   1.394 |                  1   |    nan     |     nan     | nan |      nan     |     nan     |              nan |
| window m=180: KNN vs corr RMSE |  nan     |    nan     | nan     |                  nan     | nan     |                nan   |     73.436 |      73.023 | 304 |        0.413 |       0.69  |                1 |
| window m=240: KNN vs corr RMSE |  nan     |    nan     | nan     |                  nan     | nan     |                nan   |     73.996 |      73.419 | 244 |        0.577 |       0.891 |                1 |
| window m=300: KNN vs corr RMSE |  nan     |    nan     | nan     |                  nan     | nan     |                nan   |     73.581 |      73.016 | 184 |        0.565 |       0.904 |                1 |

## G. Months 241-484 only

| W                                        | kind     |   direct |   indirect |   rho_A |   theta_hat |   se_mbb |   c_min |   equivalent |   equivalent_hac |   dyn_c_min |   dyn_equivalent |
|:-----------------------------------------|:---------|---------:|-----------:|--------:|------------:|---------:|--------:|-------------:|-----------------:|------------:|-----------------:|
| geographic KNN                           | estimate |    1.157 |      0.281 |   0.5   |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| inverse-distance                         | estimate |    1.126 |      0.58  |   0.466 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| training-correlation                     | estimate |    1.13  |      0.69  |   0.481 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| graphical-lasso                          | estimate |    1.112 |      0.342 |   0.523 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| geographic KNN vs inverse-distance       | direct   |  nan     |    nan     | nan     |       0.031 |    0.03  |   0.067 |            1 |                1 |       0.418 |                0 |
| geographic KNN vs training-correlation   | direct   |  nan     |    nan     | nan     |       0.027 |    0.032 |   0.067 |            1 |                1 |       0.412 |                0 |
| geographic KNN vs graphical-lasso        | direct   |  nan     |    nan     | nan     |       0.044 |    0.023 |   0.068 |            1 |                1 |       0.432 |                0 |
| inverse-distance vs training-correlation | direct   |  nan     |    nan     | nan     |      -0.004 |    0.009 |   0.018 |            1 |                1 |       0.159 |                0 |
| inverse-distance vs graphical-lasso      | direct   |  nan     |    nan     | nan     |       0.013 |    0.019 |   0.039 |            1 |                1 |       0.352 |                0 |
| training-correlation vs graphical-lasso  | direct   |  nan     |    nan     | nan     |       0.017 |    0.019 |   0.044 |            1 |                1 |       0.412 |                0 |
| geographic KNN vs inverse-distance       | indirect |  nan     |    nan     | nan     |      -0.299 |    0.103 |   1.662 |            0 |                0 |     nan     |              nan |
| geographic KNN vs training-correlation   | indirect |  nan     |    nan     | nan     |      -0.409 |    0.134 |   2.251 |            0 |                0 |     nan     |              nan |
| geographic KNN vs graphical-lasso        | indirect |  nan     |    nan     | nan     |      -0.061 |    0.032 |   0.505 |            0 |                0 |     nan     |              nan |
| inverse-distance vs training-correlation | indirect |  nan     |    nan     | nan     |      -0.11  |    0.048 |   0.69  |            0 |                0 |     nan     |              nan |
| inverse-distance vs graphical-lasso      | indirect |  nan     |    nan     | nan     |       0.238 |    0.095 |   1.35  |            0 |                0 |     nan     |              nan |
| training-correlation vs graphical-lasso  | indirect |  nan     |    nan     | nan     |       0.348 |    0.125 |   1.858 |            0 |                0 |     nan     |              nan |

## I. Standardised anomalies

| W                                        | kind     |   direct |   indirect |   rho_A |   phi11_mean |   theta_hat |   se_mbb |   c_min |   equivalent |   equivalent_hac |   dyn_c_min |   dyn_equivalent |
|:-----------------------------------------|:---------|---------:|-----------:|--------:|-------------:|------------:|---------:|--------:|-------------:|-----------------:|------------:|-----------------:|
| geographic KNN                           | estimate |    1.121 |      0.235 |   0.468 |        0.179 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| inverse-distance                         | estimate |    1.101 |      0.453 |   0.44  |        0.282 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| training-correlation                     | estimate |    1.116 |      0.483 |   0.424 |        0.282 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| graphical-lasso                          | estimate |    1.11  |      0.248 |   0.471 |        0.201 |     nan     |  nan     | nan     |          nan |              nan |     nan     |              nan |
| geographic KNN vs inverse-distance       | direct   |  nan     |    nan     | nan     |      nan     |       0.02  |    0.013 |   0.038 |            1 |                1 |       0.336 |                0 |
| geographic KNN vs training-correlation   | direct   |  nan     |    nan     | nan     |      nan     |       0.005 |    0.013 |   0.024 |            1 |                1 |       0.213 |                0 |
| geographic KNN vs graphical-lasso        | direct   |  nan     |    nan     | nan     |      nan     |       0.011 |    0.008 |   0.022 |            1 |                1 |       0.196 |                0 |
| inverse-distance vs training-correlation | direct   |  nan     |    nan     | nan     |      nan     |      -0.015 |    0.005 |   0.021 |            1 |                1 |       0.218 |                0 |
| inverse-distance vs graphical-lasso      | direct   |  nan     |    nan     | nan     |      nan     |      -0.009 |    0.012 |   0.026 |            1 |                1 |       0.256 |                0 |
| training-correlation vs graphical-lasso  | direct   |  nan     |    nan     | nan     |      nan     |       0.006 |    0.012 |   0.023 |            1 |                1 |       0.224 |                0 |
| geographic KNN vs inverse-distance       | indirect |  nan     |    nan     | nan     |      nan     |      -0.218 |    0.074 |   1.213 |            0 |                0 |     nan     |              nan |
| geographic KNN vs training-correlation   | indirect |  nan     |    nan     | nan     |      nan     |      -0.248 |    0.085 |   1.408 |            0 |                0 |     nan     |              nan |
| geographic KNN vs graphical-lasso        | indirect |  nan     |    nan     | nan     |      nan     |      -0.012 |    0.013 |   0.154 |            0 |                0 |     nan     |              nan |
| inverse-distance vs training-correlation | indirect |  nan     |    nan     | nan     |      nan     |      -0.03  |    0.02  |   0.257 |            0 |                0 |     nan     |              nan |
| inverse-distance vs graphical-lasso      | indirect |  nan     |    nan     | nan     |      nan     |       0.205 |    0.072 |   1.146 |            0 |                0 |     nan     |              nan |
| training-correlation vs graphical-lasso  | indirect |  nan     |    nan     | nan     |      nan     |       0.235 |    0.081 |   1.326 |            0 |                0 |     nan     |              nan |

A/H done [1s]
B done [1s]
C done [1s]
D done [1s]
E done [3s]
F done [4s]
G done [4s]
I done [4s]
