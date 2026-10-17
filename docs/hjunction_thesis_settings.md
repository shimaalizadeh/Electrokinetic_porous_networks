# H-junction: thesis settings, phd_codes input, voltage sweep, solver diff (Oct 16 2026)

## 1. Published settings (thesis Sec. 4.6, printed pp. 72-73, Fig. 4.9-4.10; Langmuir 33, 6220 (2017), Fig. 9-10)
| quantity | published | source |
|---|---|---|
| Channels 1, 2 | 6.5 mm (Lref = 6594 um in scripts) | thesis p.72; Drive hjunction/hjunct_property.py |
| Channels 3, 4 | 8.7 mm (8731 um) | thesis p.72; hjunct_property.py |
| microchannels | h = 1.4 um, w = 50 um, sigma = -13 mC/m^2 -> sigma* = -497.36, lambda* = 0.013688 | thesis p.72 |
| element | 90 um, d_n = 10 nm circular, sigma_n = -200 mC/m^2, 40% of 1 um x 50 um -> sigma* = -27.33, lambda* = 3.833, area 0.2857 | thesis p.72-73 |
| reservoirs | 1 mM, P = 0; ch.3/4 ends V = 0 | thesis p.73 |
| voltage | V_L = 10 V, V_R = 4 V (Fig. 4.9); 40 per volt in scripts -> 400 / 160 | Fig. 4.9 caption; shocks_distance.py |
| internal node volumes | largest adjacent cell (ch.1 / ch.3 dx) | thesis p.73 |
| mesh | 3200 cells per microchannel, 30 for the element | thesis p.73; paper p.6227 |
| dt | 1e-6 | thesis p.73; paper |
| run length | "long-term time advancement" to steady state; scripts use tmax = 2 | thesis p.74; binary_shock_distance.py |
| artificial diffusion | not mentioned anywhere in thesis or papers | grep of final_pdf_plain.txt, paper_part1/2.txt |
| result | d_s = 156.6 um at V_R = 4 V (d_n = 10 nm); 80.4 um for d_n = 6 nm / -140 mC/m^2 | Drive shocks_distance.py data behind Fig. 4.10 |

## 2. Is phd_codes/network_solver/input_files/hjunction/Network.in the published case? No.
- Matches: micro sigma*/lambda*, lengths 1 / 0.0136 / 1.3241, element area 0.2857, V_L = 400.
- Differs: V_R = 320 (8 V, not 4 V); element sigma* = -10.93, lambda* = 9.58 (= d 4 nm, not 10 nm);
  mesh 200/200/10/100/100 (published 3200 x4 / 30); internal volumes 0.005 / 0.0132 (published dx: 3.1e-4 / 4.1e-4);
  phd_codes main.cpp hard-codes dt = 1e-5 (main.cpp:29), published 1e-6.
- Artificial-diffusion switch (first line, read at network_info.cpp:187): 0 in Network.in, Network_other.in,
  network_solver/Network.in and single_channel/Network.in. The Drive Jan-2018 network.in (6400/80, d 6 nm) is the newer format without a switch line.
- No evidence of ramping, special initial conditions or damping in the thesis/papers: the thesis attributes stability to the implicit scheme and fine resolution (pp. 54, 67, Ch. 3).
- ArtDiff.dat in phd_codes is an empty output file.

## 3. Runs with the unchanged solver (central)
Thesis settings (examples/hjunction_reference/network.toml) RUN STABLY:
| mesh | dt | t reached | d_s (Lref 6594 um) | I_tot (last) |
|---|---|---|---|---|
| 3200/3200/30/3200/3200 | 1e-5 | 0.38 (still running) | 99.9 um, steady since t~0.3 | 400.5 |
| same | 1e-6 (thesis) | 0.035 (still running) | tracks the dt 1e-5 run (464.7 um at t=0.035) | 252.5 |
| 6400/6400/80/6400/6400 (Drive mesh) | 1e-5 | 0.186 (still running) | ~600 um, not converged | 305.8 |
Comparison: 99.9 um vs published 156.6 um (-36%). Strong mesh sensitivity (6400/80 run differs a lot), so not converged.
Figures: docs/figures/hjunction_thesis_* (see docs/figures/README.md).

The original phd_codes binary (built read-only from phd_codes, unmodified) on its own hjunction Network.in:
- switch 0, dt 1e-5: aborts ("concentration get nan") at t = 5.3e-4
- switch 1 (artificial diffusion on): aborts at t = 2.1e-4
- potentials 50% / 25% / 12%: abort at t = 5.2e-4 / 6.3e-3 / 7.4e-3
So the 2017-era solver could not run this input either; it is not a regression of the current code.

## 4. Voltage sweep, phd_codes input (both potentials scaled), current solver, dt 1e-5, tmax 0.03
See the table in the report; highest stable: 12% on the original mesh (also stable to t = 0.3, 42.7 s),
13% and above abort. The 1600/1600/10/2119/2119 mesh is stable only at <= 10%.
To reproduce: tools/voltage_sweep.py (scales both potentials of examples/hjunction_phd_codes/network.toml).
