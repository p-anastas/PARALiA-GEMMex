#!/bin/bash

# 1) Set system env variables for PARALiA-GEMMex build and execution. 
# TODO: Modify this file for your system libs
source config_system.sh
PARALIA_GEMMEX_BUILD=${PARALIA_GEMMEX_ROOT}/${PARALIA_GEMMEX_SYSTEM}-build
PARALIA_GEMMEX_DEP_FILES=${PARALIA_GEMMEX_ROOT}/Deployment_files/${PARALIA_GEMMEX_SYSTEM}

# 2) Build PARALiA-GEMMex
mkdir -p ${PARALIA_GEMMEX_BUILD} && cd ${PARALIA_GEMMEX_BUILD}
cmake ../ && make -j
cp ${PARALIA_GEMMEX_ROOT}/config_system.sh ${PARALIA_GEMMEX_DEP_FILES}/config_system_${PARALIA_GEMMEX_SYSTEM}.sh

# 3) Download and build nvbandwidth
git clone https://github.com/NVIDIA/nvbandwidth.git
mkdir -p nvbandwidth/build
cd nvbandwidth/build
cmake ../
make -j
cd ${PARALIA_GEMMEX_ROOT}
# If you want to skip step 3 and install nvbandwidth by hand (or already have it), comment out the commands above and replace NVBW_BUILD with your installation. 
NVBW_BUILD=${PARALIA_GEMMEX_BUILD}/nvbandwidth/build

# 4) Perform transfer microbenchmarks with nvbandwidth to initialize a (basic) system topology map for routing optimization
NVBW_BENCH_DIR=${PARALIA_GEMMEX_DEP_FILES}/microbench_logs
mkdir -p ${NVBW_BENCH_DIR}

${NVBW_BUILD}/nvbandwidth -t device_to_device_memcpy_write_ce > ${NVBW_BENCH_DIR}/device_to_device_memcpy_write_ce.log
${NVBW_BUILD}/nvbandwidth -t device_to_device_bidirectional_memcpy_write_ce > ${NVBW_BENCH_DIR}/device_to_device_bidirectional_memcpy_write_ce.log

${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_memcpy_ce  > ${NVBW_BENCH_DIR}/host_to_device_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_bidirectional_memcpy_ce  > ${NVBW_BENCH_DIR}/host_to_device_bidirectional_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_memcpy_ce  > ${NVBW_BENCH_DIR}/device_to_host_memcpy_ce.log
${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_bidirectional_memcpy_ce  > ${NVBW_BENCH_DIR}/device_to_host_bidirectional_memcpy_ce.log

numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_memcpy_ce  > ${NVBW_BENCH_DIR}/host_to_device_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t host_to_device_bidirectional_memcpy_ce  > ${NVBW_BENCH_DIR}/host_to_device_bidirectional_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_memcpy_ce  > ${NVBW_BENCH_DIR}/device_to_host_memcpy_ce_inter.log
numactl --interleave=all  ${NVBW_BUILD}/nvbandwidth --disableAffinity -t device_to_host_bidirectional_memcpy_ce  > ${NVBW_BENCH_DIR}/device_to_host_bidirectional_memcpy_ce_inter.log

# 5) Parse nvbandwidth logs to PARALiA-GEMMex Grid_amalgamation layout.
python3 ${PARALIA_GEMMEX_ROOT}/nvidia_topo_parse.py

# 6) Run PARALiA tests to confirm successfull installation & Grid_amalgamation layout. 
${PARALIA_GEMMEX_INSTALL_PREFIX}/testing-bin/dgemm_tester 1 1 0
${PARALIA_GEMMEX_INSTALL_PREFIX}/testing-bin/sgemm_tester 1 1 0
