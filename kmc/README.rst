

This benchmark takes two arguments on the command line. The first is the directory containing the resource files required for the simulation. The second is the input file containing the number of steps that the simulation should perform. An example launch command line can be found in the ``launch_example.sh`` script in the ``hello_world`` directory.

I suggest running this benchmark as hybrid MPI+OpenMP with 4 MPI ranks per socket.

Single Node Ballpark time estimates
-----------------------------------

Dual Gold 6142 (16 core SKL), icc 19.0.3, intel-mpi 2019 Update 3:
    8  MPI ranks 4 OMP threads: ~6.4 seconds per step

See PDF for best estimates.

