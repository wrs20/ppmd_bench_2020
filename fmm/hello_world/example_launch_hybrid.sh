#!/bin/bash
export OMP_NUM_THREADS=4
mpirun -n 4 --bind-to socket python ../fmm_strong_scaling.py ../resource_files ./hello_world.json
