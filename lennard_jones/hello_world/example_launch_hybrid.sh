#!/bin/bash
export OMP_NUM_THREADS=4
mpirun -n 4 --bind-to socket python ../lennard_jones.py ./hello_world.json
