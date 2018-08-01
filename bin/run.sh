#!/usr/bin/env bash

codedir="~/Desktop/modeling/electrokinetic_netowrk_solver"
logdir="~/Desktop/modeling/log/"
logfile="network_solver.log"
program="./main"

num_pores=774
num_reservoirs=606
hp_x=18
hp_y=33
h_pores=386
dt=1e-4
tmax=100
period=10000
restart=0
table_path="~/tables/"


mkdir -p logdir
cd $code_dir

###########################RUN PROGRAM###############################
$program ${num_pores} ${num_reservoirs} ${hp_x} ${hp_y}  ${h_pores} ${dt} ${tmax}  ${period}  ${restart} ${table_path} &> ${logdir}/logfile

