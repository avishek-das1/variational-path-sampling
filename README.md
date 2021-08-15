# Direct evaluation of rare events in active matter from variational path sampling
This repository contains the basic code for VPS and the data for reproducing the corresponding paper [arXiv:2108.05359](arxiv.org/abs/2108.05359).

The codes directory contains optimization routines for well-tempered metadynamics, MCVB-T and MCVB algorithms. Other than the files already supplied, it requires a collection of 10000 initial conditions for trajectories, to be called 'v9_frames.txt'. We have only supplied one frame as an example due to file-size constraints. Compilation and execution commands:

mpifort -O3 mcvb.f90 -o mcvb

mpirun -np 2 ./mcvb

The figure_data directory contains the data and corresponding python scripts to plot the figures. For generating 'fig1.png', command is 

python fig1.py
