# H-junction from phd_codes (network_solver/input_files/hjunction)

- `network.toml`: converted from the original v0 `Network.in` + `Channels.in`:
  `python3 tools/network_in_to_toml.py Network.in -p 5 -r 6 --v0 --channels Channels.in`
  (200/200/10/100/100 cells). Blows up in the current solver within ~4e-5 (dt 1e-5 to 1e-7).
- `network_fine.toml`: same physics, mesh refined to 6400/6400/10/8475/8475 cells
  (dx ~1.56e-4 in the micro-channels; the 10-cell nano-bridge unchanged). This runs:

      mkdir -p run/input_files && cp examples/hjunction_phd_codes/network_fine.toml run/input_files/
      ln -s ../tables run/tables && cd run
      ../build/ek_solver 5 6 5 1e-5 0.03 30 0 tables/ network_fine.toml > log.txt   # ~3.5 min
      cd .. && ekpn-viz hjunction run --network run/input_files/network_fine.toml \
          --dt 1e-5 --write-every 30 --log run/log.txt --out figs
- `network_recommended.toml`: result of the mesh study. It is the fine mesh. With the unchanged
  solver, stability is not monotonic in mesh refinement, so there is no coarser mesh that runs reliably.
- Voltage sweep: with potentials scaled to 12% (48 / 38.4) the original mesh runs stably; 13% and above abort
  (tools/voltage_sweep.py). NOTE: this input is not the published case (V_R = 8 V, 4 nm element); see
  examples/hjunction_reference/ for the published settings.
