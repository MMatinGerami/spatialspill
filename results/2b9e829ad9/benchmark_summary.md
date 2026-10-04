# Simulator benchmark (geometry: perturb_fish)

| scenario        | estimator   | kind       |   ring |   null_fpr_q |   null_cov95 |   power_q |   auroc |   fdp_q |   n_tests |
|:----------------|:------------|:-----------|-------:|-------------:|-------------:|----------:|--------:|--------:|----------:|
| auto_only       | E1          | autonomous |     -1 |            0 |        0.913 |         0 |   0.936 |     nan |       900 |
| auto_only       | E1          | spillover  |      1 |            0 |        0.908 |       nan | nan     |     nan |       150 |
| auto_only       | E1          | spillover  |      2 |            0 |        0.912 |       nan | nan     |     nan |       900 |
| everything      | E1          | spillover  |      2 |            0 |        0.898 |         0 |   0.496 |     nan |       330 |
| no_effect       | E1          | autonomous |     -1 |            0 |        0.958 |       nan | nan     |     nan |       900 |
| no_effect       | E1          | spillover  |      1 |            0 |        0.893 |       nan | nan     |     nan |       150 |
| no_effect       | E1          | spillover  |      2 |            0 |        0.918 |       nan | nan     |     nan |       900 |
| spill           | E1          | autonomous |     -1 |            0 |        0.89  |         0 |   0.913 |     nan |       900 |
| spill           | E1          | spillover  |      1 |            0 |        0.879 |         0 |   0.399 |     nan |       150 |
| spill           | E1          | spillover  |      2 |            0 |        0.909 |         0 |   0.499 |     nan |       900 |
| spill_batch     | E1          | autonomous |     -1 |            0 |        0.897 |         0 |   0.902 |     nan |       900 |
| spill_batch     | E1          | spillover  |      1 |            0 |        0.899 |         0 |   0.459 |     nan |       150 |
| spill_batch     | E1          | spillover  |      2 |            0 |        0.911 |         0 |   0.515 |     nan |       900 |
| spill_bleed     | E1          | autonomous |     -1 |            0 |        0.886 |         0 |   0.913 |     nan |       900 |
| spill_bleed     | E1          | spillover  |      1 |            0 |        0.846 |         0 |   0.444 |     nan |       150 |
| spill_bleed     | E1          | spillover  |      2 |            0 |        0.924 |         0 |   0.504 |     nan |       900 |
| spill_clonal    | E1          | spillover  |      2 |            0 |        0.886 |         0 |   0.38  |     nan |       270 |
| spill_density   | E1          | autonomous |     -1 |            0 |        0.885 |         0 |   0.905 |     nan |       900 |
| spill_density   | E1          | spillover  |      1 |            0 |        0.896 |         0 |   0.313 |     nan |       150 |
| spill_density   | E1          | spillover  |      2 |            0 |        0.915 |         0 |   0.422 |     nan |       900 |
| spill_misassign | E1          | autonomous |     -1 |            0 |        0.904 |         0 |   0.912 |     nan |       900 |
| spill_misassign | E1          | spillover  |      1 |            0 |        0.88  |         0 |   0.612 |     nan |       210 |
| spill_misassign | E1          | spillover  |      2 |            0 |        0.916 |         0 |   0.532 |     nan |       900 |
