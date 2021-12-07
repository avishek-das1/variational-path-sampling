# Direct evaluation of rare events in active matter from variational path sampling
This repository contains the basic code for VPS and the data for reproducing the corresponding paper [arXiv:2108.05359](arxiv.org/abs/2108.05359).

The codes directory contains routines for generating uncorrelated initial snapshots and then performing well-tempered metadynamics, MCVB-T and MCVB optimization algorithms. The initial snapshots file must first be generated as the other routines need it. Compilation and execution commands are as follows.

For generating initial snapshots:

mpifort -O3 initialize.f90

./a.out

For optimization:

mpifort -O3 mcvb.f90 -o mcvb

mpirun -np 2 ./mcvb

The figure_data directory contains the data and corresponding python scripts to plot the figures. For generating 'fig1.png', command is 

python fig1.py
