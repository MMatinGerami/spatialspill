# Simulator benchmark (geometry: perturb_fish)

| scenario        | estimator   | kind       |   ring |   null_fpr_q |   null_cov95 |   power_q |   auroc |   fdp_q |   n_tests |
|:----------------|:------------|:-----------|-------:|-------------:|-------------:|----------:|--------:|--------:|----------:|
| auto_only       | E1          | autonomous |     -1 |        0.019 |        0.911 |     0.759 |   0.936 |   0.176 |       900 |
| auto_only       | E1          | spillover  |      2 |        0.01  |        0.951 |   nan     | nan     |   1     |       900 |
| everything      | E1          | spillover  |      2 |        0.023 |        0.93  |     0     |   0.526 |   1     |       330 |
| no_effect       | E1          | autonomous |     -1 |        0.001 |        0.961 |   nan     | nan     |   1     |       900 |
| no_effect       | E1          | spillover  |      2 |        0.004 |        0.947 |   nan     | nan     |   1     |       900 |
| spill           | E1          | autonomous |     -1 |        0.026 |        0.891 |     0.726 |   0.913 |   0.237 |       900 |
| spill           | E1          | spillover  |      2 |        0.015 |        0.943 |     0     |   0.521 |   1     |       900 |
| spill_batch     | E1          | autonomous |     -1 |        0.025 |        0.902 |     0.684 |   0.901 |   0.239 |       900 |
| spill_batch     | E1          | spillover  |      2 |        0.016 |        0.944 |     0.083 |   0.518 |   0.972 |       900 |
| spill_bleed     | E1          | autonomous |     -1 |        0.028 |        0.891 |     0.758 |   0.913 |   0.24  |       900 |
| spill_bleed     | E1          | spillover  |      2 |        0.014 |        0.948 |     0     |   0.514 |   1     |       900 |
| spill_clonal    | E1          | spillover  |      2 |        0.042 |        0.907 |     0     |   0.394 |   1     |       270 |
| spill_density   | E1          | autonomous |     -1 |        0.024 |        0.886 |     0.748 |   0.903 |   0.219 |       900 |
| spill_density   | E1          | spillover  |      2 |        0.016 |        0.948 |     0     |   0.431 |   1     |       900 |
| spill_misassign | E1          | autonomous |     -1 |        0.019 |        0.905 |     0.716 |   0.912 |   0.188 |       900 |
| spill_misassign | E1          | spillover  |      2 |        0.016 |        0.938 |     0     |   0.547 |   1     |       900 |
