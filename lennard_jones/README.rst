
This benchmark takes one argument on the command line. This argument is the input file containing the number of steps that the simulation should perform. An example launch command line can be found in the ``example_launch_hybrid.sh`` script in the ``hello_world`` directory.

I suggest running this benchmark as hybrid MPI+OpenMP with 4 MPI ranks per socket.

Single Node Ballpark time estimates
-----------------------------------

Dual Gold 6142 (16 core SKL), icc 19.0.3, intel-mpi 2019 Update 3:
    8  MPI ranks 4 OMP threads: ~0.042 seconds per step
