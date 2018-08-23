#!/usr/bin/env bash

codedir="/Users/alizadeh/Desktop/Electrokinetic_porous_networks"
logdir="/Users/alizadeh/Desktop/Electrokinetic_porous_networks"
logfile="network_solver.log"
program="./main"

num_pores=1
num_reservoirs=2
h_pores=1
dt=1e-4
tmax=1
period=1000000
restart=0
table_path="tables/"
network_input_file="Network.in"

cd ${codedir}

###########################RUN PROGRAM###############################
$program ${num_pores} ${num_reservoirs} ${h_pores} ${dt} ${tmax} ${period} ${restart} ${table_path} ${network_input_file}

