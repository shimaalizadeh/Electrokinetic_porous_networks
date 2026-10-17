# Mesh study, unchanged central solver (tools/mesh_study.py)

Stable = reaches max_time without abort and I_tot within 0.5% of the finest mesh (sp: 3200 cells, tj: 12800 uniform, hj: 6400/6400/10/8475/8475 fine run, 121.335).
min C_bar is taken over ~200 snapshots, so it is approximate. Wall time is on 8 shared cores with up to 8 runs in parallel.

| case | mesh (cells per pore) | dt | result | t_end | I_tot | dev vs finest | min C_bar | wall [s] |
|---|---|---|---|---|---|---|---|---|
| H-junction phd_codes | 400/400/10/530/530 | 1e-05 | ABORT | 4e-05 | aborted | - | 0.000902 | 0.0 |
| H-junction phd_codes | 800/800/10/1059/1059 | 1e-05 | ABORT | 7e-05 | aborted | - | 0.000902 | 0.1 |
| H-junction phd_codes | 1600/1600/10/2119/2119 | 1e-05 | ABORT | 0.0002 | aborted | - | -1.07 | 0.4 |
| H-junction phd_codes | 3200/3200/10/4237/4237 | 1e-05 | ABORT | 0.00018 | aborted | - | 0.000902 | 0.7 |
| H-junction phd_codes | 4800/4800/10/6356/6356 | 1e-05 | ABORT | 0.00022 | aborted | - | 0.000902 | 1.4 |
| H-junction phd_codes | 5600/5600/10/7415/7415 | 1e-05 | ABORT | 0.00017 | aborted | - | -0.453 | 1.2 |
| H-junction phd_codes | 6000/6000/10/7945/7945 | 1e-05 | ABORT | 0.00016 | aborted | - | -2.19 | 1.4 |
| H-junction phd_codes | 6200/6200/10/8209/8209 | 1e-05 | ok | 0.03001 | 121.3360 | +0.00% | -0.973 | 186.9 |
| H-junction phd_codes | 6400/6400/5/8474/8474 | 1e-05 | ABORT | 0.0003 | aborted | - | -2.66 | 2.3 |
| H-junction phd_codes | 6400/6400/10/8474/8474 | 1e-05 | ok | 0.03001 | 121.3350 | +0.00% | -0.676 | 186.8 |
| H-junction phd_codes | 6400/6400/20/8474/8474 | 1e-05 | ABORT | 0.00029 | aborted | - | -0.000239 | 2.2 |
| H-junction phd_codes | 6600/6600/10/8739/8739 | 1e-05 | ABORT | 0.00019 | aborted | - | -0.456 | 1.8 |
| H-junction phd_codes | 6800/6800/10/9004/9004 | 1e-05 | ok | 0.03001 | 121.3350 | +0.00% | -0.287 | 205.3 |
| H-junction phd_codes | 8000/8000/10/10593/10593 | 1e-05 | ABORT | 0.00031 | aborted | - | -6.85 | 2.9 |
| H-junction phd_codes | 9600/9600/10/12711/12711 | 1e-05 | ABORT | 0.00018 | aborted | - | -0.719 | 2.2 |
| H-junction phd_codes | 12800/12800/10/16948/16948 | 1e-05 | ABORT | 0.00031 | aborted | - | -0.418 | 4.6 |
| H-junction phd_codes | 800/800/10/1059/1059 | 1e-06 | ABORT | 4.1e-05 | aborted | - | 0.000902 | 0.4 |
| H-junction phd_codes | 3200/3200/10/4237/4237 | 2e-06 | ABORT | 0.000214 | aborted | - | 0.000902 | 3.7 |
| H-junction phd_codes | 6400/6400/10/8474/8474 | 5e-06 | ABORT | 0.00015 | aborted | - | -430 | 2.2 |
| H-junction phd_codes | 6400/6400/10/8474/8474 | 2e-05 | ABORT | 0.00014 | aborted | - | -0.574 | 0.8 |
| single pore (V=80) | 100 | 0.0001 | ABORT | 0.0002 | aborted | - | 0.975 | 0.1 |
| single pore (V=80) | 200 | 0.0001 | ABORT | 0.0006 | aborted | - | -3.06 | 0.1 |
| single pore (V=80) | 300 | 0.0001 | ABORT | 0.0006 | aborted | - | -32 | 0.1 |
| single pore (V=80) | 350 | 0.0001 | ABORT | 0.0008 | aborted | - | -0.0437 | 0.1 |
| single pore (V=80) | 400 | 0.0001 | ok | 0.1 | 23.2686 | -0.23% | -0.0235 | 1.2 |
| single pore (V=80) | 800 | 0.0001 | ok | 0.1 | 23.3003 | -0.10% | -0.00192 | 1.9 |
| single pore (V=80) | 1600 | 0.0001 | ok | 0.1 | 23.3155 | -0.03% | -0.00118 | 3.4 |
| single pore (V=80) | 3200 | 0.0001 | ok | 0.1 | 23.3230 | +0.00% | -2e-05 | 6.7 |
| T-junction V=80 | 800/800/800 | 0.0001 | ABORT | 0.0005 | aborted | - | 0.0251 | 0.1 |
| T-junction V=80 | 1600/1600/1600 | 0.0001 | ABORT | 0.0006 | aborted | - | -2.54 | 0.2 |
| T-junction V=80 | 800/2400/2400 | 0.0001 | ok | 0.05 | 55.7791 | +0.07% | -0.00456 | 7.2 |
| T-junction V=80 | 2000/2000/2000 | 0.0001 | ABORT | 0.0008 | aborted | - | -44.8 | 0.2 |
| T-junction V=80 | 1200/2400/2400 | 0.0001 | ok | 0.05 | 55.7439 | +0.01% | -0.00695 | 7.8 |
| T-junction V=80 | 3200/1600/1600 | 0.0001 | ABORT | 0.0007 | aborted | - | -1.08 | 0.2 |
| T-junction V=80 | 1600/2400/2400 | 0.0001 | ok | 0.05 | 55.7340 | -0.01% | -0.00729 | 8.6 |
| T-junction V=80 | 2400/2400/2400 | 0.0001 | ok | 0.05 | 55.7291 | -0.02% | -0.00734 | 8.4 |
| T-junction V=80 | 1600/3200/3200 | 0.0001 | ok | 0.05 | 55.7360 | -0.00% | -0.00381 | 9.4 |
| T-junction V=80 | 2800/2800/2800 | 0.0001 | ok | 0.05 | 55.7298 | -0.02% | -0.00533 | 9.9 |
| T-junction V=80 | 3200/3200/3200 | 0.0001 | ok | 0.05 | 55.7306 | -0.01% | -0.00421 | 11.9 |
| T-junction V=80 | 6400/6400/6400 | 0.0001 | ok | 0.05 | 55.7351 | -0.01% | -0.00372 | 23.3 |
| T-junction V=80 | 12800/12800/12800 | 0.0001 | ok | 0.05 | 55.7384 | +0.00% | -0.0038 | 45.9 |
| T-junction V=80 | 1600/1600/1600 | 1e-05 | ABORT | 0.00034 | aborted | - | 0.0505 | 0.4 |
| T-junction V=80 | 6400/6400/6400 | 1e-05 | ok | 0.05 | 55.7351 | -0.01% | -0.0045 | 214.2 |

## Recommendations
- Single pore: 400 cells (examples/single_pore/network_recommended.toml). 350 and below abort.
- T-junction at V=80: 800/2400/2400 (inlet pore / junction branches), dt 1e-4
  (examples/tjunction/network_v80_recommended.toml). The coarsest uniform mesh that runs is 2400. The branches
  that meet at the junction set the limit; refining only the inlet pore (3200/1600/1600) still aborts.
- H-junction (phd_codes): 6400/6400/10/8475/8475 at dt 1e-5, which is network_fine.toml
  (examples/hjunction_phd_codes/network_recommended.toml). Stability is not monotonic in the mesh:
  f = 31, 32 and 34 ran, while f = 30, 33, 40, 48 and 64 aborted, as did smaller and larger dt and bridge 5 or 20.
  So there is no robust "coarsest stable" mesh. The working meshes agree on I_tot (121.335-121.336), but
  the bridge C_bar goes to -0.3 to -1.0, so the runs that survive are not clean either.

