#!/bin/bash

## This is the configuration template file for PARALiA-GEMMex compilation and deployment in a new system.
## The intended use is to create one such config for each different deplyed system (e.g. config_sys1.sh, config_sys2.sh etc) in order to allow easy deployment.
#------------------------------------------------------Required Modules/bindings----------------------------------------------------#
module load CUDA/12.6.0 #NVHPC-HPCX-CUDA12/24.9
module load Boost/1.82.0-GCC-12.3.0
module load OpenBLAS/0.3.23-GCC-12.3.0
module load Python/3.11.3-GCCcore-12.3.0
module load CMake/3.26.3-GCCcore-12.3.0
source ~/.bashrc
export OMP_PROC_BIND=spread
#--------------------------------------------------------------Basic----------------------------------------------------------------#

# CHECKME: A desired name for the testbed to be used for your build-dirs and logfiles.
system="vega_4A100"
export PARALIA_SYSTEM=${system}

# CHECKME: A folder to change to before firing benchmarks, to avoid disk errors. 
export PARALIA_CD_FOLDER=/ceph/hpc/home/eupetrosa/

# Define the (max) num of devices PARALiA can use (will pick the first 'num' dev_ids)
export PARALIA_NUM_DEVICES=4

# CHECKME: Define cuda architecture (Tesla K40 = 35, GTX 1060/70 = 61) P100 = 60, V100 = 70, A100 = 80
export PARALIA_CUDA_ARCH=80

# CHECKME: Define cuda path. Use "default" to use the default compiler calls without prefix(es)
export PARALIA_CUDA_TOOLKIT_PREFIX="default"

# CHECKME: Define cuda load command (e.g. should always be '-lcuda', but its not! If you don't know what to do, leave it)
export PARALIA_CUDA_LOAD_COMMAND='-lcuda'

# CHECKME: Define gcc compiler path. Use "default" to use the default compiler calls without prefix(es)
export PARALIA_CXX_PREFIX="default"

# CHECKME: Define path for prebuild openblas. NOTE: OpenBLAS built using the same gcc is adviced.
export PARALIA_OPENBLAS_PREFIX="default"

# CHECKME: Define path for prebuild boost. NOTE: boost built using the same gcc is adviced.
export PARALIA_BOOST_PREFIX="default"

# CHECKME: Define the directory the project will be installed in. 'default' = buildir/${system}-install
export PARALIA_INSTALL_PREFIX="default"

# CHECKME: Define the Watt for your CPU - no CPU power tool used currently.
export PARALIA_W_CPU_PREDEF=560

#-----------------------------------------------------------------------------------------------------------------------------------#
