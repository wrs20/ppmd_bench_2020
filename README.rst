Installation
------------

In this directory there should be a ``INSTALL.rst`` that describes how to setup an environment for running these benchmarks.


Benchmarks
----------

1. Lennard-Jones: This benchmark represents a model Molecular Dynamics system, it typically reachs a strong scaling limit around 100 nodes.
2. FMM: A benchmark that implements a Fast Multipole method. This is an example of a electrostatic solver that uses a hierarchical grid structure.
3. KMC: A Kinetic Monte Carlo benchmark. This benchmark should scale very well as the ratio of computation to communication is high. This benchmark should be compute bound.

Output
------

Each run of a benchmark code should be performed in a separate directory. On completion of a run the codes will produce a json output file ``last_time.json`` that contains the keys:

``time_taken_per_step``
    The time taken per step of simulation, this is our observable of interest.

``num_omp_threads``
    Recorded for data collection quality of life.

``num_mpi_ranks``
    Recorded for data collection quality of life.


