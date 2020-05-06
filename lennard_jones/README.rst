
I suggest running this benchmark as hybrid MPI+OpenMP with 4 MPI ranks per socket.

Single Node Ballpark time estimates
-----------------------------------

Dual Gold 6142 (16 core SKL), icc 19.0.3, intel-mpi 2019 Update 3:
    8  MPI ranks 4 OMP threads: ~0.042 seconds per step
