Introduction
------------

This benchmark takes two inputs on the command line. The first is the directory containing the resource files required for the simulation. The second is the input file containing the number of steps that the simulation should perform. An example launch command line can be found in the "hello_world.sh" script in the "hello_world" directory.

This benchmark has a hard strong scaling limit of 4096 MPI ranks and will not launch with more. Typically I place 1 or 2 MPI ranks per socket and use hybrid MPI+OpenMP.


Ballpark time estimates
-----------------------

Dual E5-2680, clang 7.0.0, open-mpi 1.10.2:
    4  MPI rank  4 OMP threads: ~11.5 seconds per step
    16 MPI ranks 1 OMP thread : ~10.7 seconds per step 


