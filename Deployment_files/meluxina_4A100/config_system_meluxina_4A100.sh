#!/bin/bash

## This is the configuration template file for PARALiA-GEMMex compilation and deployment in a new system.
## The intended use is to create one such config for each different deplyed system (e.g. config_sys1.sh, config_sys2.sh etc) in order to allow easy deployment.

## NOTE: This file must be source'd every time for PARALiA-GEMMex correctly!!!

#------------------------------------------------------Required Modules/bindings----------------------------------------------------#
# TODO: Here you have to load all required modules for PARALiA-GEMMex cmake/compilation/execution
# NOTE: Skip this step if working in a system that is not module-based (might require including your own paths/installations in the next step) 

# 1. A CUDA compiler (Tested 10.X,11.X,12.X) or an SDK that contains one
module load CUDA
module load numactl

# 2. Boost
module load Boost/1.82.0-GCC-12.3.0

# 3. OpenBLAS
module load OpenBLAS/0.3.23-GCC-12.3.0

# 4. Cmake
module load CMake/3.26.3-GCCcore-12.3.0

# 5. Python (this is not required for running the library, but it is for plotters and deployment helper scripts)
module load Python/3.11.3-GCCcore-12.3.0

# 6. Process/OMP/GPU binding commands or configurations for efficient execution on each HPC system.
source ~/.bashrc
export OMP_PROC_BIND=spread
#--------------------------------------------------------------Basic----------------------------------------------------------------#

# TODO: A desired name for the testbed to be used for your build-dirs and logfiles.
export PARALIA_GEMMEX_SYSTEM="meluxina_4A100"

# TODO: Change this if you run deploy from another dir
export PARALIA_GEMMEX_ROOT=$(pwd)

# TODO: Define the directory the project will be installed in. 'default' = ${PARALIA_GEMMEX_SYSTEM}-build/${PARALIA_GEMMEX_SYSTEM}-install
export PARALIA_GEMMEX_INSTALL_PREFIX=${PARALIA_GEMMEX_ROOT}/${PARALIA_GEMMEX_SYSTEM}-build/${PARALIA_GEMMEX_SYSTEM}-install

# TODO: Define cuda architecture (Tesla K40 = 35, GTX 1060/70 = 61) P100 = 60, V100 = 70, A100 = 80
export PARALIA_GEMMEX_CUDA_ARCH=80

# TODO: Define the (max) num of devices PARALiA-GEMMex can use (will pick the first 'num' dev_ids). 
# Usually = the number of GPUs in the system.
export PARALIA_GEMMEX_NUM_DEVICES=4

#--------------------------------------------------------------Extra----------------------------------------------------------------#

# CHECKME: Define cuda path. Leave "default" to (try to) use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CUDA_PREFIX="default"

# CHECKME: Define cuda load command (e.g. should always be '-lcuda', but sometimes its not! If you don't know what to do, leave it)
export PARALIA_GEMMEX_CUDA_LOAD_COMMAND='-lcuda'

# CHECKME: Define gcc compiler path. Use "default" to (try to) use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CXX_PREFIX="default"

# CHECKME: Define path for prebuild openblas. NOTE: OpenBLAS built using the same gcc is adviced. Most modules cover this, but not all...
export PARALIA_GEMMEX_OPENBLAS_PREFIX="default"

# CHECKME: Define path for prebuild boost. NOTE: boost built using the same gcc is adviced. Most modules cover this, but not all...
export PARALIA_GEMMEX_BOOST_PREFIX="default"

#-----------------------------------------------------------------------------------------------------------------------------------#
