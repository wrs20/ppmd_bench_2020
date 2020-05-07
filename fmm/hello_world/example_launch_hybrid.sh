#!/bin/bash

# For 2 nodes with 2 x 8 core cpus cores 
export OMP_NUM_THREADS=2
mpirun -ppn 8 -np 16  python ../fmm_strong_scaling.py ../resource_files ./hello_world.json
