#!/bin/bash

## This is the configuration template file for PARALiA-GEMMex compilation and deployment in a new system.
## The intended use is to create one such config for each different deplyed system (e.g. config_sys1.sh, config_sys2.sh etc) in order to allow easy deployment.

#------------------------------------------------------Required Modules/bindings----------------------------------------------------#
export OMP_NUM_THREADS=20
export OMP_PROC_BIND=spread
#--------------------------------------------------------------Basic----------------------------------------------------------------#

# CHECKME: A desired name for the testbed to be used for your build-dirs and logfiles.
system="silver1_2V100"
export PARALIA_GEMMEX_SYSTEM=${system}

# CHECKME: A folder to change to before firing benchmarks, to avoid disk errors. 
export PARALIA_GEMMEX_CD_FOLDER=/mnt/proj1/dd-23-129

# Define the (max) num of devices PARALiA can use (will pick the first 'num' dev_ids)
export PARALIA_GEMMEX_NUM_DEVICES=2

# CHECKME: Define cuda architecture (Tesla K40 = 35, GTX 1060/70 = 61,) P100 = 60, V100 = 70, A100 = 80
export PARALIA_GEMMEX_CUDA_ARCH=70

# CHECKME: Define cuda toolkit path. If left 'default' PARALiA-GEMMex will try to use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CUDA_PREFIX="/usr/local/cuda-12.1"

# CHECKME: Define cuda load command (e.g. should always be '-lcuda', but its not! If you don't know what to do, leave it)
export PARALIA_GEMMEX_CUDA_LOAD_COMMAND='-lcuda'

# CHECKME: Define gcc compiler path. If left 'default' PARALiA-GEMMex will try to use the default compiler calls without prefix(es)
export PARALIA_GEMMEX_CXX_PREFIX="default"

# CHECKME: Define path for prebuild openblas. NOTE: OpenBLAS built using the same gcc is adviced.
export PARALIA_GEMMEX_OPENBLAS_PREFIX="/home/users/panastas/Lib_install/OpenBLAS-gcc93"

# CHECKME: Define path for prebuild boost. NOTE: boost built using the same gcc is adviced.
export PARALIA_GEMMEX_BOOST_PREFIX="/home/users/panastas/Lib_install/boost-1.69"

# CHECKME: Define the directory PARALiA-GEMMex will be installed in. 'default' = buildir/${system}-install
export PARALIA_GEMMEX_INSTALL_PREFIX="default"
#-----------------------------------------------------------------------------------------------------------------------------------#
