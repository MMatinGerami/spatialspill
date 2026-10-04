# Simulator benchmark (geometry: perturb_fish)

| scenario        | estimator   | kind       |   ring |   null_fpr_q |   null_cov95 |   power_q |   auroc |   fdp_q |   n_tests |
|:----------------|:------------|:-----------|-------:|-------------:|-------------:|----------:|--------:|--------:|----------:|
| auto_only       | E1          | autonomous |     -1 |        0.02  |        0.913 |     0.732 |   0.936 |   0.189 |       900 |
| auto_only       | E1          | spillover  |      1 |        0.006 |        0.908 |   nan     | nan     |   1     |       150 |
| auto_only       | E1          | spillover  |      2 |        0.009 |        0.912 |   nan     | nan     |   1     |       900 |
| everything      | E1          | spillover  |      2 |        0.009 |        0.898 |     0     |   0.496 |   1     |       330 |
| no_effect       | E1          | autonomous |     -1 |        0     |        0.958 |   nan     | nan     | nan     |       900 |
| no_effect       | E1          | spillover  |      1 |        0     |        0.893 |   nan     | nan     | nan     |       150 |
| no_effect       | E1          | spillover  |      2 |        0     |        0.918 |   nan     | nan     | nan     |       900 |
| spill           | E1          | autonomous |     -1 |        0.019 |        0.89  |     0.726 |   0.913 |   0.182 |       900 |
| spill           | E1          | spillover  |      1 |        0.013 |        0.879 |     0     |   0.399 |   1     |       150 |
| spill           | E1          | spillover  |      2 |        0.005 |        0.909 |     0     |   0.499 |   1     |       900 |
| spill_batch     | E1          | autonomous |     -1 |        0.017 |        0.897 |     0.689 |   0.902 |   0.174 |       900 |
| spill_batch     | E1          | spillover  |      1 |        0.028 |        0.899 |     0     |   0.459 |   1     |       150 |
| spill_batch     | E1          | spillover  |      2 |        0.003 |        0.911 |     0     |   0.515 |   1     |       900 |
| spill_bleed     | E1          | autonomous |     -1 |        0.023 |        0.886 |     0.726 |   0.913 |   0.214 |       900 |
| spill_bleed     | E1          | spillover  |      1 |        0.01  |        0.846 |     0     |   0.444 |   1     |       150 |
| spill_bleed     | E1          | spillover  |      2 |        0.003 |        0.924 |     0     |   0.504 |   1     |       900 |
| spill_clonal    | E1          | spillover  |      2 |        0.004 |        0.886 |     0     |   0.38  |   1     |       270 |
| spill_density   | E1          | autonomous |     -1 |        0.019 |        0.885 |     0.732 |   0.905 |   0.182 |       900 |
| spill_density   | E1          | spillover  |      1 |        0.01  |        0.896 |     0     |   0.313 |   1     |       150 |
| spill_density   | E1          | spillover  |      2 |        0.006 |        0.915 |     0     |   0.422 |   1     |       900 |
| spill_misassign | E1          | autonomous |     -1 |        0.03  |        0.904 |     0.732 |   0.912 |   0.26  |       900 |
| spill_misassign | E1          | spillover  |      1 |        0.013 |        0.88  |     0     |   0.612 |   1     |       210 |
| spill_misassign | E1          | spillover  |      2 |        0.012 |        0.916 |     0     |   0.532 |   1     |       900 |
