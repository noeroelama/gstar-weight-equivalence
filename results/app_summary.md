# Application results (generated 2026-10-04 04:48)

Panel 484 x 40; scored 38; window m=240; B=1999; alpha=0.05; relative margin 0.1 of geographic KNN; delta_RMSE=1.0 mm

## In-sample estimates (anomalies, GSTARX)

| W                    |   phi10_mean |   phi11_mean |   phi11_min |   phi11_max |   n_phi11_negative |   direct |   indirect |   total |   rho_A |   ljung_box_reject_share |   mean_abs_resid_crosscorr |   resid_sd_mm |   indirect_h3 |   indirect_h6 |   indirect_h12 |   indirect_h24 |   se_phi11_mean |   se_direct |   se_indirect |
|:---------------------|-------------:|-------------:|------------:|------------:|-------------------:|---------:|-----------:|--------:|--------:|-------------------------:|---------------------------:|--------------:|--------------:|--------------:|---------------:|---------------:|----------------:|------------:|--------------:|
| geographic KNN       |        0.095 |        0.161 |      -0.246 |       0.561 |                  5 |    1.155 |      0.214 |   1.369 |   0.456 |                    0.375 |                      0.273 |        66.701 |         0.211 |         0.215 |          0.214 |          0.214 |           0.029 |       0.023 |         0.05  |
| inverse-distance     |        0.096 |        0.243 |      -0.119 |       0.578 |                  3 |    1.135 |      0.397 |   1.532 |   0.401 |                    0.35  |                      0.276 |        66.621 |         0.375 |         0.396 |          0.397 |          0.397 |           0.049 |       0.021 |         0.108 |
| training-correlation |        0.1   |        0.265 |      -0.061 |       0.557 |                  3 |    1.134 |      0.474 |   1.608 |   0.401 |                    0.275 |                      0.276 |        66.494 |         0.437 |         0.471 |          0.474 |          0.474 |           0.049 |       0.02  |         0.127 |
| graphical-lasso      |        0.068 |        0.196 |      -0.297 |       0.722 |                  6 |    1.123 |      0.257 |   1.38  |   0.487 |                    0.325 |                      0.274 |        66.673 |         0.253 |         0.258 |          0.257 |          0.257 |           0.032 |       0.019 |         0.057 |

## Equivalence Q1/Q2

| functional   | W_A                  | W_B                  |   theta_hat |   se_hac |   se_mbb |   L_pw |   L_used |   ci90_lo |   ci90_hi |   z_se |   min_delta_abs |   tau_ref |   delta_ref |   c_min_ref | equivalent_ref   |   c_min_ref_a025 | equivalent_ref_a025   |   c_min_ref_hac | equivalent_ref_hac   |   tau_dyn |   delta_dyn |   c_min_dyn |   equivalent_dyn |   c_min_dyn_a025 |   equivalent_dyn_a025 |   c_min_dyn_hac |   equivalent_dyn_hac |
|:-------------|:---------------------|:---------------------|------------:|---------:|---------:|-------:|---------:|----------:|----------:|-------:|----------------:|----------:|------------:|------------:|:-----------------|-----------------:|:----------------------|----------------:|:---------------------|----------:|------------:|------------:|-----------------:|-----------------:|----------------------:|----------------:|---------------------:|
| phi11_mean   | geographic KNN       | inverse-distance     |     -0.0821 |   0.0282 |   0.0264 | 2.787  |        3 |   -0.1256 |   -0.0386 | 0.0435 |          0.1256 |    0.1612 |      0.0161 |      0.7365 | False            |           0.7845 | False                 |          0.7395 | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| phi11_mean   | geographic KNN       | training-correlation |     -0.1037 |   0.0304 |   0.0285 | 2.6372 |        3 |   -0.1506 |   -0.0569 | 0.0468 |          0.1506 |    0.1612 |      0.0161 |      0.941  | False            |           1.01   | False                 |          0.929  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| phi11_mean   | geographic KNN       | graphical-lasso      |     -0.0344 |   0.011  |   0.0103 | 2.6372 |        3 |   -0.0514 |   -0.0174 | 0.017  |          0.0514 |    0.1612 |      0.0161 |      0.3425 | False            |           0.374  | False                 |          0.3475 | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| phi11_mean   | inverse-distance     | training-correlation |     -0.0216 |   0.0091 |   0.0091 | 2.6372 |        3 |   -0.0367 |   -0.0066 | 0.015  |          0.0367 |    0.1612 |      0.0161 |      0.2725 | False            |           0.3085 | False                 |          0.263  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| phi11_mean   | inverse-distance     | graphical-lasso      |      0.0477 |   0.0254 |   0.0248 | 3.3175 |        4 |    0.0068 |    0.0885 | 0.0409 |          0.0885 |    0.1612 |      0.0161 |      0.503  | False            |           0.5425 | False                 |          0.505  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| phi11_mean   | training-correlation | graphical-lasso      |      0.0693 |   0.0263 |   0.0246 | 2.8774 |        3 |    0.0288 |    0.1098 | 0.0405 |          0.1098 |    0.1612 |      0.0161 |      0.678  | False            |           0.7335 | False                 |          0.667  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| direct       | geographic KNN       | inverse-distance     |      0.02   |   0.0154 |   0.0143 | 2.7478 |        3 |   -0.0036 |    0.0435 | 0.0235 |          0.0435 |    1.1549 |      0.1155 |      0.0375 | True             |           0.041  | True                  |          0.039  | True                 |    0.1549 |      0.0155 |      0.2675 |                0 |           0.2945 |                     0 |          0.2765 |                    0 |
| direct       | geographic KNN       | training-correlation |      0.0205 |   0.0165 |   0.0153 | 2.7478 |        3 |   -0.0047 |    0.0456 | 0.0252 |          0.0456 |    1.1549 |      0.1155 |      0.039  | True             |           0.043  | True                  |          0.0405 | True                 |    0.1549 |      0.0155 |      0.2765 |                0 |           0.303  |                     0 |          0.281  |                    0 |
| direct       | geographic KNN       | graphical-lasso      |      0.0318 |   0.0113 |   0.011  | 2.7478 |        3 |    0.0138 |    0.0499 | 0.0181 |          0.0499 |    1.1549 |      0.1155 |      0.0425 | True             |           0.0455 | True                  |          0.043  | True                 |    0.1549 |      0.0155 |      0.302  |                0 |           0.321  |                     0 |          0.301  |                    0 |
| direct       | inverse-distance     | training-correlation |      0.0005 |   0.0058 |   0.0058 | 2.7478 |        3 |   -0.009  |    0.01   | 0.0095 |          0.01   |    1.1549 |      0.1155 |      0.009  | True             |           0.0105 | True                  |          0.009  | True                 |    0.1549 |      0.0155 |      0.0675 |                1 |           0.0835 |                     1 |          0.0675 |                    1 |
| direct       | inverse-distance     | graphical-lasso      |      0.0119 |   0.0136 |   0.0128 | 2.9639 |        3 |   -0.0091 |    0.0328 | 0.021  |          0.0328 |    1.1549 |      0.1155 |      0.0285 | True             |           0.032  | True                  |          0.03   | True                 |    0.1549 |      0.0155 |      0.2215 |                0 |           0.2525 |                     0 |          0.233  |                    0 |
| direct       | training-correlation | graphical-lasso      |      0.0113 |   0.0132 |   0.0127 | 2.7478 |        3 |   -0.0096 |    0.0323 | 0.021  |          0.0323 |    1.1549 |      0.1155 |      0.0285 | True             |           0.032  | True                  |          0.029  | True                 |    0.1549 |      0.0155 |      0.223  |                0 |           0.2565 |                     0 |          0.2355 |                    0 |
| indirect     | geographic KNN       | inverse-distance     |     -0.1831 |   0.0703 |   0.0667 | 3.5355 |        4 |   -0.2928 |   -0.0734 | 0.1097 |          0.2928 |    0.2143 |      0.0214 |      1.2175 | False            |           1.297  | False                 |          1.236  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| indirect     | geographic KNN       | training-correlation |     -0.2596 |   0.0908 |   0.0865 | 3.1483 |        4 |   -0.4018 |   -0.1174 | 0.1422 |          0.4018 |    0.2143 |      0.0214 |      1.712  | False            |           1.8285 | False                 |          1.7165 | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| indirect     | geographic KNN       | graphical-lasso      |     -0.0431 |   0.0178 |   0.0173 | 1.3878 |        2 |   -0.0716 |   -0.0146 | 0.0285 |          0.0716 |    0.2143 |      0.0214 |      0.3835 | False            |           0.435  | False                 |          0.373  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| indirect     | inverse-distance     | training-correlation |     -0.0765 |   0.0263 |   0.0237 | 1.3878 |        2 |   -0.1155 |   -0.0375 | 0.039  |          0.1155 |    0.2143 |      0.0214 |      0.5745 | False            |           0.634  | False                 |          0.55   | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| indirect     | inverse-distance     | graphical-lasso      |      0.14   |   0.0648 |   0.0606 | 3.7496 |        4 |    0.0403 |    0.2397 | 0.0997 |          0.2397 |    0.2143 |      0.0214 |      0.981  | False            |           1.048  | False                 |          0.999  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |
| indirect     | training-correlation | graphical-lasso      |      0.2165 |   0.0831 |   0.0785 | 3.3773 |        4 |    0.0874 |    0.3456 | 0.1291 |          0.3456 |    0.2143 |      0.0214 |      1.4315 | False            |           1.5225 | False                 |          1.448  | False                |  nan      |    nan      |    nan      |              nan |         nan      |                   nan |        nan      |                  nan |

## RMSE

| W                    | model   |   rmse_mm |   M |
|:---------------------|:--------|----------:|----:|
| geographic KNN       | GSTAR   |    75.014 | 244 |
| inverse-distance     | GSTAR   |    74.082 | 244 |
| training-correlation | GSTAR   |    73.716 | 244 |
| graphical-lasso      | GSTAR   |    74.99  | 244 |
| geographic KNN       | GSTARX  |    73.996 | 244 |
| inverse-distance     | GSTARX  |    73.687 | 244 |
| training-correlation | GSTARX  |    73.419 | 244 |
| graphical-lasso      | GSTARX  |    73.996 | 244 |

## Equivalence RMSE

| model   | W_A                  | W_B                  |   theta_rmse_mm |   delta_mm |   L_pw |   L_used |   se_mbb |   se_hac | equivalent_a025   |   min_delta_mm_a025 |   ci90_lo |   ci90_hi | equivalent_mbb   | equivalent_hac   |   p_tost_mbb |   min_delta_mm |   se_mbb_L2 |   se_mbb_L6 |   se_mbb_L12 |   se_mbb_L24 |
|:--------|:---------------------|:---------------------|----------------:|-----------:|-------:|---------:|---------:|---------:|:------------------|--------------------:|----------:|----------:|:-----------------|:-----------------|-------------:|---------------:|------------:|------------:|-------------:|-------------:|
| GSTAR   | geographic KNN       | inverse-distance     |           0.932 |          1 |  2.953 |        3 |    0.225 |    0.24  | False             |               1.372 |     0.562 |     1.301 | False            | False            |        0.38  |          1.301 |       0.218 |       0.244 |        0.2   |        0.199 |
| GSTAR   | geographic KNN       | training-correlation |           1.298 |          1 |  3.39  |        4 |    0.342 |    0.341 | False             |               1.969 |     0.735 |     1.861 | False            | False            |        0.808 |          1.861 |       0.307 |       0.347 |        0.314 |        0.299 |
| GSTAR   | geographic KNN       | graphical-lasso      |           0.024 |          1 |  2.227 |        3 |    0.076 |    0.074 | True              |               0.173 |    -0.102 |     0.149 | True             | True             |        0     |          0.149 |       0.075 |       0.073 |        0.076 |        0.08  |
| GSTAR   | inverse-distance     | training-correlation |           0.366 |          1 |  3.593 |        4 |    0.125 |    0.126 | True              |               0.611 |     0.161 |     0.571 | True             | True             |        0     |          0.571 |       0.109 |       0.128 |        0.132 |        0.139 |
| GSTAR   | inverse-distance     | graphical-lasso      |          -0.908 |          1 |  3.088 |        4 |    0.237 |    0.24  | False             |               1.373 |    -1.299 |    -0.518 | False            | False            |        0.349 |          1.299 |       0.217 |       0.243 |        0.209 |        0.197 |
| GSTAR   | training-correlation | graphical-lasso      |          -1.274 |          1 |  3.553 |        4 |    0.329 |    0.336 | False             |               1.92  |    -1.816 |    -0.732 | False            | False            |        0.797 |          1.816 |       0.293 |       0.349 |        0.304 |        0.306 |
| GSTARX  | geographic KNN       | inverse-distance     |           0.309 |          1 |  1.873 |        2 |    0.121 |    0.123 | True              |               0.547 |     0.11  |     0.509 | True             | True             |        0     |          0.509 |       0.123 |       0.119 |        0.114 |        0.108 |
| GSTARX  | geographic KNN       | training-correlation |           0.577 |          1 |  2.72  |        3 |    0.188 |    0.196 | True              |               0.946 |     0.267 |     0.887 | True             | True             |        0.012 |          0.887 |       0.185 |       0.2   |        0.203 |        0.209 |
| GSTARX  | geographic KNN       | graphical-lasso      |          -0     |          1 |  1     |        2 |    0.062 |    0.064 | True              |               0.122 |    -0.102 |     0.101 | True             | True             |        0     |          0.102 |       0.06  |       0.064 |        0.06  |        0.062 |
| GSTARX  | inverse-distance     | training-correlation |           0.268 |          1 |  3.479 |        4 |    0.089 |    0.092 | True              |               0.442 |     0.122 |     0.414 | True             | True             |        0     |          0.414 |       0.084 |       0.092 |        0.103 |        0.112 |
| GSTARX  | inverse-distance     | graphical-lasso      |          -0.31  |          1 |  1     |        2 |    0.11  |    0.105 | True              |               0.524 |    -0.49  |    -0.13  | True             | True             |        0     |          0.49  |       0.112 |       0.103 |        0.092 |        0.08  |
| GSTARX  | training-correlation | graphical-lasso      |          -0.578 |          1 |  2.106 |        3 |    0.172 |    0.171 | True              |               0.915 |    -0.861 |    -0.295 | True             | True             |        0.007 |          0.861 |       0.165 |       0.175 |        0.176 |        0.176 |

## Data facts

|                        |       0 |
|:-----------------------|--------:|
| mean_rain_mm           | 224.873 |
| mean_rain_scored_mm    | 222.098 |
| sd_ratio_median        |   2.415 |
| sd_ratio_min           |   1.585 |
| sd_ratio_max           |  33.558 |
| sd_ratio_median_scored |   2.415 |

## Raw data (no seasonal adjustment)

| W                    |   phi10_mean |   phi11_mean |   direct |   indirect |   rho_A |   ljung_box_reject_share |
|:---------------------|-------------:|-------------:|---------:|-----------:|--------:|-------------------------:|
| geographic KNN       |        0.286 |        0.196 |    1.738 |      0.248 |   0.871 |                        1 |
| inverse-distance     |        0.247 |        0.37  |    1.45  |      1.241 |   0.667 |                        1 |
| training-correlation |        0.312 |        0.289 |    1.571 |      0.998 |   0.671 |                        1 |
| graphical-lasso      |        0.235 |        0.268 |    1.62  |      0.407 |   0.835 |                        1 |

in-sample fits done [2s]
Q1/Q2 equivalence done [2s]
rolling forecasts done [4s]
RMSE equivalence done [4s]
