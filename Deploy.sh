#!/bin/bash

# 1) Set system env variables for PARALiA-GEMMex build and execution. 
# TODO: Modify this file for your system libs
source config_system.sh

# 2) Build PARALiA-GEMMex
mkdir ${PARALIA_GEMMEX_SYSTEM}_build
cd ${PARALIA_GEMMEX_SYSTEM}_build
cmake ../
make -j

# 3) Download and build nvbandwidth
git clone https://github.com/NVIDIA/nvbandwidth.git
mkdir -p nvbandwidth/build
cd nvbandwidth/build
cmake ../
make -j
cd ../../../
# If you want to skip step 3 and install nvbandwidth by hand (or already have it), comment out the commands above and replace NVBW_BUILD with your installation. 
NVBW_BUILD=${PARALIA_GEMMEX_SYSTEM}_build/nvbandwidth/build

# 4) Download and build nvbandwidth
mkdir -p Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs

${NVBW_BUILD}/nvbandwidth -t device_to_device_memcpy_write_ce > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_device_memcpy_write_ce.log
${NVBW_BUILD}/nvbandwidth -t device_to_device_bidirectional_memcpy_write_ce > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_device_bidirectional_memcpy_write_ce.log

${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/host_to_device_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_bidirectional_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/host_to_device_bidirectional_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_host_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_bidirectional_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_host_bidirectional_memcpy_ce.log

numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/host_to_device_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_bidirectional_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/host_to_device_bidirectional_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_host_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_bidirectional_memcpy_ce  > Deployment_files/${PARALIA_GEMMEX_SYSTEM}/microbench_logs/device_to_host_bidirectional_memcpy_ce_inter.log

python3 ../nvidia_topo_parse.py