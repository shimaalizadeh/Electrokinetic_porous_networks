# H-junction with the published thesis settings (Sec. 4.6, Langmuir 2017 Part II)
- network.toml: 3200 cells per microchannel, 30 for the element (thesis); run with dt 1e-6 (thesis) or 1e-5 (same d_s).
- network_6400.toml: 6400/80 mesh of the Jan-2018 Drive network.in (but with the thesis 10 nm element).
- Shock distance: python3 tools/hjunction_ds.py RUN --dt DT --write-every N. Details: docs/hjunction_thesis_settings.md.
