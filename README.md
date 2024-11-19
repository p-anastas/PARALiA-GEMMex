# PARALiA-GEMMex (Uncut-GEMMs)

This repository provides a high-performance General Matrix Multiplication (GEMM) library for single-node multi-GPU HPC systems. The library uses autotuning to analyze each problem's communication, scheduling, and imbalance characteristics on the fly and create a customized static schedule with optimized data movement, caching, and overlap for that problem (more in [our paper](https://www.computer.org/csdl/proceedings-article/cluster/2024/587100a143/21HPrlbPFzG)).

## Software/Compiler requirements
 - CUDA toolkit 10+ (Latest release tested with 12.x versions on V100 and A100 clusters)
 - A gcc/g++ compiler compatible with the above CUDA (tested with 11.x, 12.x)
 - OpenBLAS, installed with the same gcc compiler.
 - Boost, installed with the same gcc compiler.
 - CMake minimum version 3.10
 - numactl
 - Python 3.x (packages: os, pandas. Additionally for plotting: math, numpy, scipy, matplotlib, seaborn)
 - The [nvbandwidth](https://github.com/NVIDIA/nvbandwidth) tool (Installation & microbenchmarks already included in `./deploy.sh`, see bellow).
## Installation
PARALiA-GEMMex installation consists of 2(3) easy steps:
 - Fill `config_system.sh` with compiler/library details, module loads etc.
 - (Optional) Modify `CMakeLists.txt` to enable/disable optimizations and/or check experimental features.
 - Run `./deploy.sh` *on the target system* (cross-compilation not supported due to deployment microbenchmarks). This file:
   - Installs PARALiA-GEMMex.
   - Downloads and installs nvbandwidth (if done manually, you have to set its path instead).
   - Performs some fast microbenchmarks (~mins) to estimate the communication capabilities of the target interconnect.
   - Runs tests for PARALiA-GEMMex functions to validate installation correctness (~mins).

## Usage
After a succesful installation, you should have:
  - `${PARALIA_GEMMEX_INSTALL_PREFIX}/lib`, that contains shared .so files (Main library functions: `libparalia.so`)
  - `${PARALIA_GEMMEX_INSTALL_PREFIX}/include`,  that contains all header files (Main header functions: `PARALiA.hpp`)
  - `${PARALIA_GEMMEX_INSTALL_PREFIX}/testing-bin` with tests and benchmarks for _PARALiA-GEMMex_ and _cuBLASXt_ routines. 

To use PARALiA-GEMMex functions:
 - Always: `source config_system.sh` 
 - Compiling: Include the main header `PARALiA.hpp` in code and  `-I${PARALIA_GEMMEX_INSTALL_PREFIX}/include` during compilation.
 - Linking: use `-L${PARALIA_GEMMEX_INSTALL_PREFIX}/lib -lparalia` during linking.
 - PARALiA BLAS functions accept usual BLAS paramaters in a similar way to OpenBLAS/cuBLASXt etc (also see benchmarks bellow).
   - `PARALiADgemm(TransA, TransB, M, N, K, alpha, A, ldA, B, ldB, beta, C, ldC)`
     - _A, B, C_ can reside in CPU or (any) GPU memory.

## Prebuild Benchmarks
 - `${PARALIA_GEMMEX_INSTALL_PREFIX}/testing-bin` contains `PARALiA-GEMMex` runners for _double_ and _single_ precision
   - Also _half_, but still experimental (validators not updated for low precision).
 - Usage: `./[s,d]gemm_runner dev_num dev_ids T cache_max_size TransA TransB alpha beta D1 D2 D3 loc1 loc2 loc3 outloc`
   - **Control parameters** (you should leave all these to _-1_ unless you know what you are doing):
     - `dev_num`: The number of GPUs for the benchmark. Use all system GPUs if < 0.
     - `dev_ids`: A list of the GPUs used for execution. Ignored if dev_num < 0. 
       - Input form example: `0101 for devices = [0,2]`, `1111 for devices = [0,1,2,3]` etc.
     - `T`: The internal tiling size. Problem-tailored automatic tiling size selection if < 0.
     - `cache_max_size`: The maximum cache size that can be allocated in each GPU. Defined automatically if < 0.
   - **Routine input parameters**:
     - `TransA, TransB`: **_N_** or **_T_**. The transpose parameter used for GEMM routine invocation for the A,B matrices, respectively.
     - `alpha, beta`: The GEMM constants for the routine invocation.
     - `D1 D2 D3`: D1 = M, D2 = N, D3 = K for the routine invocation.
    - **Data placement parameters**:
     - `loc1 loc2 loc3`: The locations of the A,B and C matrices in memory. Input form example:
       - `loc = 0 to (system_gpu_num - 1)`: Matrix initially on the corresponding **GPU memory** (order = nvidia-smi).
       - `loc = system_gpu_num`: Matrix on pinned **host memory** (Numa-unaware, not advised).
       - `loc = system_gpu_num + 1`: Matrix on pinned **numa-interleaved host memory** (advised for host allocations).
     - `outloc`: The output location of the C matrix. **Always** set to `loc3` currently.

## Contact
For questions, issues, or collaboration inquiries, please reach out via **email:** [panastas@cslab.ece.ntua.gr](mailto:panastas@cslab.ece.ntua.gr)

## Licence 
- This work is open-source and distributed under a GPL-3.0 license.

## Related Publications/Citations:
 - Main PARALiA-GEMMex implementation: [Uncut-GEMMs: Communication-Aware Matrix Multiplication on Multi-GPU Nodes](https://www.computer.org/csdl/proceedings-article/cluster/2024/587100a143/21HPrlbPFzG)
 - Scheduling & autotuning: [PARALiA: A Performance Aware Runtime for Auto-tuning Linear Algebra on heterogeneous systems](https://dl.acm.org/doi/10.1145/3624569)
 - Overlap modeling & T selection: [CoCoPeLia: Communication-Computation Overlap Prediction for Efficient Linear Algebra on GPUs](https://ieeexplore.ieee.org/document/9408195)

