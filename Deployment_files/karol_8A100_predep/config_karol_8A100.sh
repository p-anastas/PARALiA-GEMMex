#!/bin/bash

## This is the configuration template file for PARALiA-GEMMex compilation and deployment in a new system.
## The intended use is to create one such config for each different deplyed system (e.g. config_sys1.sh, config_sys2.sh etc) in order to allow easy deployment.
#------------------------------------------------------Required Modules/bindings----------------------------------------------------#
module load Boost/1.81.0-GCC-12.2.0
module load OpenBLAS/0.3.21-GCC-12.2.0
module load Python/3.10.8-GCCcore-12.2.0
module load CMake/3.24.3-GCCcore-12.2.0
module load NVHPC/24.1-CUDA-12.4.0

#source ~/.bashrc
export OMP_PROC_BIND=spread
#--------------------------------------------------------------Basic----------------------------------------------------------------#

# CHECKME: A desired name for the testbed to be used for your build-dirs and logfiles.
system="karol_8A100"
export PARALIA_GEMMEX_SYSTEM=${system}

# CHECKME: A folder to change to before firing benchmarks, to avoid disk errors. 
export PARALIA_GEMMEX_CD_FOLDER=/mnt/proj1/dd-23-129

# Define the (max) num of devices PARALiA can use (will pick the first 'num' dev_ids)
export PARALIA_GEMMEX_NUM_DEVICES=8

# CHECKME: Define cuda architecture (Tesla K40 = 35, GTX 1060/70 = 61,) P100 = 60, V100 = 70, A100 = 80
export PARALIA_GEMMEX_CUDA_ARCH=80

# CHECKME: Define cuda toolkit path. Use "default" to use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CUDA_PREFIX="/usr/local/cuda-12.4"

# CHECKME: Define cuda load command (e.g. should always be '-lcuda', but its not! If you don't know what to do, leave it)
export PARALIA_GEMMEX_CUDA_LOAD_COMMAND='-lcuda'

# CHECKME: Define gcc compiler path. Use "default" to use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CXX_PREFIX="default"

# CHECKME: Define path for prebuild openblas. NOTE: OpenBLAS built using the same gcc is adviced.
export PARALIA_GEMMEX_OPENBLAS_PREFIX="/apps/all/OpenBLAS/0.3.21-GCC-12.2.0"

# CHECKME: Define path for prebuild boost. NOTE: boost built using the same gcc is adviced.
export PARALIA_GEMMEX_BOOST_PREFIX="/apps/all/Boost/1.81.0-GCC-12.2.0"

# CHECKME: Define the directory the project will be installed in. 'default' = buildir/${system}-install
export PARALIA_GEMMEX_INSTALL_PREFIX="default"

#-----------------------------------------------------------------------------------------------------------------------------------#
