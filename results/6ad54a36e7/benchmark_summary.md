# Simulator benchmark (geometry: synthetic)

| scenario   | estimator    | kind       |   ring |   null_fpr_q |   null_cov95 |   ntc_cov95 |   ntc_fpr_p05 |   power_q |   auroc |   fdp_q |   n_tests |
|:-----------|:-------------|:-----------|-------:|-------------:|-------------:|------------:|--------------:|----------:|--------:|--------:|----------:|
| no_effect  | E0_neighbour | autonomous |     -1 |        0.001 |        0.936 |       0.938 |         0.045 |   nan     | nan     |   1     |      1500 |
| no_effect  | E0_neighbour | spillover  |      0 |        0.011 |        0.902 |       0.828 |         0.111 |   nan     | nan     |   1     |      1080 |
| no_effect  | E2           | autonomous |     -1 |        0.002 |        0.93  |       0.937 |         0.063 |   nan     | nan     |   1     |      1200 |
| no_effect  | E2           | spillover  |      2 |        0.017 |        0.85  |     nan     |       nan     |   nan     | nan     |   1     |       120 |
| spill      | E0_neighbour | autonomous |     -1 |        0.011 |        0.94  |       0.943 |         0.038 |     0.588 |   0.894 |   0.13  |      1500 |
| spill      | E0_neighbour | spillover  |      0 |        0.01  |        0.89  |       0.861 |         0.067 |     0     |   0.447 |   1     |      1080 |
| spill      | E2           | autonomous |     -1 |        0.011 |        0.942 |       0.93  |         0.07  |     0.368 |   0.795 |   0.125 |      1200 |
| spill      | E2           | spillover  |      2 |        0.043 |        0.889 |     nan     |       nan     |     0     |   0.761 |   1     |       120 |
