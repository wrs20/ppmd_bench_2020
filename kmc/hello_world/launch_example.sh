#!/bin/bash
export OMP_NUM_THREADS=1
mpirun -n 16 --bind-to socket python ../kmc_strong_scaling.py ../resource_files hello_world.json
