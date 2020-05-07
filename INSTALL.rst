Introduction
~~~~~~~~~~~~

Instructions that set up a PPMD environment on Balena using intel compiler/mpi.


Installation
~~~~~~~~~~~~

* A recent Python 3 release (3.4 - 3.8). I am using 3.8.1

::

    module load intel/compiler/64/19.5.281
    module load intel/mpi/64/18.4.274

    # My python build
    module load ~/scratch/apps/python-3.8.1/modulefile


To install a fixed benchmark version of ``PPMD``/``coulomb_kmc`` into a Python virtual environment:

::

    pip install --upgrade --no-cache-dir git+https://github.com/ppmd/ppmd@kmc_bench_2019
    pip install --upgrade --no-cache-dir git+https://github.com/ppmd/coulomb_kmc@bench_2019
    
    
Runtime
~~~~~~~
   
These are the environment variables I currently define.

::
    
    # stops each python processes hitting the FS
    export PYTHONDONTWRITEBYTECODE=1

    # ensures that the hashes of python objects are the same on each process
    export PYTHONHASHSEED=123

    # define where PPMD should place compiled binaries. Must be reachable by all processes.
    # This directory needs to be different per architecture
    export PPMD_BUILD_DIR=/beegfs/scratch/user/m/wrs20/venvs/$VIRTUAL_ENV/build

    # Define which compilers to use. The default compilers are located in ``ppmd/config/compilers``
    # The environment variable should be set to a value which matches the ``name`` field of the
    # desired compiler. Additional compilers can be defined by the user, see docs for more info.
    
    # For intel platforms where icc xHost is reliable
    export PPMD_CC_MAIN=ICC
    export PPMC_CC_OMP=ICC
    
    # For AMD Rome plaforms
    export PPMD_CC_MAIN=ICC_ROME
    export PPMC_CC_OMP=ICC_ROME

    # disable cuda explictly (disables attempted module loading on platforms without cude devices)
    export PPMD_DISABLE_CUDA=1





    
