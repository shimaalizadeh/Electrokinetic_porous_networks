# Figures

All figures here were made with `ekpn-viz` (see the top-level README for installation) from runs of the
unchanged solver. Each section lists the commands used. Total size is kept below ~2 MB on purpose;
regenerate anything else locally.

## Single pore: deionization shock (`single_pore_*`)
Case: `input_files/Network.in` unchanged (one slit pore, L = 0.6, 400 cells, sigma* = -200, lambda_D = 0.017;
inlet reservoir at potential 80 with fixed P and mu; outlet is a small closed reservoir, volume 8.068e-5).
Same physics as the I_tot = 23.2686 baseline, with frequent snapshots over the first t = 0.05.

    mkdir -p run/input_files && cp input_files/Network.in run/input_files/ && ln -s ../tables run/tables
    cd run && ../build/ek_solver 1 2 1 1e-4 0.05 2 0 tables/ Network.in > log.txt && cd ..
    ekpn-viz single-pore run --network run/input_files/Network.in --dt 1e-4 --write-every 2 \
        --log run/log.txt --out docs/figures --mp4

- `single_pore_schematic.png`: drawing of a charged pore between two reservoirs
  (`python -m ekpn_viz.schematic docs/figures/single_pore_schematic.png`)
- `single_pore_spacetime.png`: C_bar(x, t), x across and t up, black = depleted, dashed line = front
  (C_bar = C_bar_inlet / 2); layout follows Alizadeh, PhD dissertation, Stanford University, 2017
- `single_pore_profiles.png`: profiles at log-spaced times
- `single_pore_shock.mp4` / `.gif`: the front, animated

## H-junction (`hjunction_*`)
Run: `examples/hjunction_reference/network.toml` (3200/3200/30/3200/3200 cells), dt 1e-5, V_L = 400, V_R = 160,
bottom channel grounded. Snapshots every 50 steps to t = 0.03 plus every 1000 steps to t = 0.75 (two runs merged
by symlinking Tecplot files under one write_every = 50 numbering).

    R=run; N=examples/hjunction_reference/network.toml; C="--network $N --dt 1e-5 --write-every 50 --log $R/log.txt"
    ekpn-viz animate $R $C --prefix hjunction_slit --paths "0,1;3,4" --shape slit --title "H-junction, slit channels"
    ekpn-viz render3d $R --network $N --dt 1e-5 --write-every 50 --times 0.75 --prefix hjunction_3d

- `hjunction_slit_C_bar.mp4` / `.gif`: the network drawn as slit channels coloured by C_bar (channel lengths to scale, bridge drawn short), animated, with C_bar along both channels
- `hjunction_3d_C_bar_t0.75.png`: render3d (PyVista tubes, log colour scale) at t = 0.75

## T-junction (`tjunction_*`)
`examples/tjunction/network.toml`: three 0.6-long slit pores, 400 cells, sigma* = -200, meeting at internal
node 3; V = 20 at node 0, 0 at nodes 1 and 2; dt 1e-4, t 0.05, write_every 5.

    ekpn-viz network run --network examples/tjunction/network.toml --dt 1e-4 --write-every 5 --prefix tjunction

- `tjunction_C_bar_map.png`: network coloured by C_bar
