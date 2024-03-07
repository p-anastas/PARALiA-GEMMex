#!/bin/bash

## This is the configuration template file for CHANGELINK compilation and deployment in a new system.
## The intended use is to create one such config for each different deplyed system (e.g. config_sys1.sh, config_sys2.sh etc) in order to allow easy deployment.

#--------------------------------------------------------------Basic----------------------------------------------------------------#

# CHECKME: A desired name for the testbed to be used for your build-dirs and logfiles.
system="silver1_2V100"
export PARALIA_SYSTEM=${system}

# CHECKME: A folder to change to before firing benchmarks, to avoid disk errors. 
export PARALIA_CD_FOLDER=/mnt/proj1/dd-23-129

# CHECKME: A file that contains information found by numactl and nvidia-smi in a tidied form.
# FIXME: Currently must be defined by the user, will be automated later
export PARALIA_SCONFIG_FILE_PATH="/home/users/panastas/PhD_stuff/PARALiA-GEMMex/sconfig_silver1.in"

# CHECKME: Define cuda architecture (Tesla K40 = 35, GTX 1060/70 = 61,) P100 = 60, V100 = 70, A100 = 80
export PARALIA_CUDA_ARCH=70

# CHECKME: Define cuda toolkit path. If left 'default' CHANGELINK will try to use the default compiler calls without prefix(es)
export PARALIA_CUDA_TOOLKIT_PREFIX="/usr/local/cuda-12.1"

# CHECKME: Define cuda load command (e.g. should always be '-lcuda', but its not! If you don't know what to do, leave it)
export PARALIA_CUDA_LOAD_COMMAND='-lcuda'

# CHECKME: Define gcc compiler path. If left 'default' CHANGELINK will try to use the default compiler calls without prefix(es)
export PARALIA_CXX_PREFIX="default"

# CHECKME: Define path for prebuild openblas. NOTE: OpenBLAS built using the same gcc is adviced.
export PARALIA_OPENBLAS_PREFIX="/home/users/panastas/Lib_install/OpenBLAS-gcc93"

# CHECKME: Define path for prebuild boost. NOTE: boost built using the same gcc is adviced.
export PARALIA_BOOST_PREFIX="/home/users/panastas/Lib_install/boost-1.69"

# CHECKME: Define the directory CHANGELINK will be installed in. 'default' = buildir/${system}-install
export PARALIA_INSTALL_PREFIX="default"

# CHECKME: Define the Watt for your CPU - no CPU power tool used currently.
export PARALIA_W_CPU_PREDEF=280

#-----------------------------------------------------------------------------------------------------------------------------------#