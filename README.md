# Multi-Scale Modeling of Electrokinetic Porous Networks
The developed pipeline allows efficient modeling of electrokinetic phenomena in porous networks.
A porous network is modeled as a random network of micro-scale and nano-scale pores.

# Running Code
The modeling package is implemented in C++. You need C++ compiler (such as gcc) to compile the code and build the executable file.
Make sure that you have already installed lapack and blas library on your laptop/PC.

Once you have all the requirements ready, compile codes by running the Makefile:
make clean  
make

If everything goes through, it means all things are installed properly and now you are ready to run the model!

# Input Files
To run the code you have to provide Network.in located at ./input_files folder. Please refer to the guideline on how to prepare the input files at ./input_files folder.
The developed model uses a number of look-up tables to complete computations. These tables are already generated and located at the default folder ./tables.
If you put them i na different location, you just need to specify the path in the bash script in ./bin folder. 
Once the input files are ready, a bash script is used to run the code, and it's located at bin folder.

# Toy Problems
A few example problems were provided along with source code. You can find the corresponding input files 
under ./examples folder. 

[//]: # (the input file is designed for a network that involves a loop, as shown below:)

![alt text](https://github.com/shimaalizadeh/Electrokinetic_porous_networks/blob/master/toy_problem.png)



## Network file: network.toml (new) or Network.in (legacy)

The run config's `[io] network_file` picks the reader by extension: `.toml` uses the
labeled format, anything else the legacy column format. Both fill the same solver
arrays, so results are identical (checked on 5 example networks, see `examples/`).

- Annotated schema: `input_files/network_schema.toml`
- Example: `input_files/network.toml` + `run_single_pore_network_toml.toml` (I_tot = 23.2686)
- Convert old files: `python3 tools/network_in_to_toml.py input_files/Network.in -p <num_pores> -r <num_reservoirs> -o input_files/network.toml`
- Per-cell properties (old Channels.in): add `profile = "cells.csv"` to a pore; CSV header
  `sigma_star,lambda_D,Sp,C0`, one row per cell.
