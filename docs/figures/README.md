# Figures

All figures here were made with `ekpn-viz` (see the top-level README for installation) from runs of the
unchanged solver. Each section lists the commands used. Total size is kept small (about 7 MB) on purpose;
regenerate anything else locally.
All case parameters below (lengths, potentials, times, sigma*, lambda_D, volumes) are in the solver's dimensionless
units; see "Units" in the top-level README.

## Single pore: deionization shock (`single_pore_*`)
Case: `input_files/Network.in` unchanged (one slit pore, L = 0.6, 400 cells, sigma* = -200, lambda_D = 0.017;
inlet reservoir at potential 80 with fixed P and mu; the other end is a dead end: a small closed reservoir,
volume 8.068e-5, with fixed mu and free pressure. It acts as an ideal cation-selective membrane, so the anion
flux and the net flow there are zero and the pressure builds up to oppose the electroosmotic flow).
Same physics as the I_tot = 23.2686 baseline, with frequent snapshots over the first t = 0.05.

    mkdir -p run/input_files && cp input_files/Network.in run/input_files/ && ln -s ../tables run/tables
    cd run && ../build/ek_solver 1 2 1 1e-4 0.05 2 0 tables/ Network.in > log.txt && cd ..
    ekpn-viz single-pore run --network run/input_files/Network.in --dt 1e-4 --write-every 2 \
        --log run/log.txt --out docs/figures --mp4 --stop-at-stall

`--stop-at-stall` ends the single-pore space-time plot and the video when the front first comes within 1% of the pore
length of its final position (here x = 0.055 at t = 0.032); the stalled period after that is left out.

- `single_pore_schematic.png`: drawing of the charged pore, inlet reservoir on the left and an ideal cation-selective
  membrane at the dead end on the right (zero net flow, zero anion flux), as in the dissertation Fig. 4.1a
  (`python -m ekpn_viz.schematic docs/figures/single_pore_schematic.png`)
- `single_pore_spacetime.png`: C_bar(x, t), x across and t up, black = depleted, dashed line = front
  (C_bar = C_bar_inlet / 2); layout follows Alizadeh, PhD dissertation, Stanford University, 2017
- `single_pore_thickness.png`: shock thickness vs t, where thickness = C_bar_inlet / |dC_bar/dx| at the half
  level (the tangent there, run from C_bar = 0 to C_bar_inlet), with a t^1/2 fit; curve and fit both end at the peak thickness (t = 0.0206). The curve starts once the shock has
  formed (depleted side below 10% of the inlet value, from t = 0.0004); in the first snapshot after the uniform
  initial condition (t = 0.0002) depletion is still building up at the dead end and there is no shock yet
- `single_pore_shock.mp4` / `.gif`: the front, animated

## H-junction (`hjunction_*`)
Run: `examples/hjunction_reference/network.toml` (3200/3200/30/3200/3200 cells), dt 1e-5, V_L = 400, V_R = 160,
bottom channel grounded. The slit animation uses `examples/hjunction_reference/run.toml` (snapshots every 50 steps
to t = 0.1). The 3D render used snapshots every 50 steps to t = 0.03 plus every 1000 steps to t = 0.75 (two runs
merged by symlinking Tecplot files under one write_every = 50 numbering).

    R=run; N=examples/hjunction_reference/network.toml; C="--network $N --dt 1e-5 --write-every 50 --log $R/log.txt"
    ekpn-viz animate $R $C --prefix hjunction_slit --paths "0,1;3,4" --shape slit --style clean --title "H-junction" \
        --path-names "Top channel (V = 400 to 160);Bottom channel (grounded ends)" \
        --hatch-pores 2 --pore-label "2:Cation-selective element"
    ekpn-viz render3d $R --network $N --dt 1e-5 --write-every 50 --times 0.75 --prefix hjunction_3d

- `hjunction_slit_C_bar.mp4` / `.gif`: the two channels and the cation-selective element drawn as slit channels
  colored by C_bar (channel lengths to scale, element drawn short), animated, with C_bar along both channels below
- `hjunction_3d_C_bar_t0.75.png`: render3d (PyVista tubes, log color scale) at t = 0.75

## T-junction (`tjunction_*`)
`examples/tjunction/network.toml`: three 0.6-long slit pores, 400 cells, sigma* = -200, meeting at internal
node 3; V = 20 at node 0, 0 at nodes 1 and 2; dt 1e-4, t 0.05, write_every 5.

    ekpn-viz network run --network examples/tjunction/network.toml --dt 1e-4 --write-every 5 --prefix tjunction

- `tjunction_C_bar_map.png`: network colored by C_bar
