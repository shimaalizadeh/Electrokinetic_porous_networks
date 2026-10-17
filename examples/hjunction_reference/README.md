# H-junction with the published settings (Sec. 4.6, Langmuir 2017 Part II)
- network.toml: 3200 cells per microchannel, 30 for the element (as published); run with dt 1e-6 (as published) or 1e-5 (same d_s).
- network_6400.toml: 6400/80 mesh of the Jan-2018 Drive network.in (but with the published 10 nm element).
- Shock distance: python3 tools/hjunction_ds.py RUN --dt DT --write-every N. Details: docs/hjunction_thesis_settings.md.
