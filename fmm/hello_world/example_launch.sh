#!/bin/bash
export OMP_NUM_THREADS=1
mpirun -n 16 --bind-to none python ../fmm_strong_scaling.py ../resource_files ./hello_world.json
