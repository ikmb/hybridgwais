# HybridGWAIS

## Requirements
- C++ compiler with C++17 support
- Boost library
- pthread
- OpenMP
- zlib
- Intel TBB >= 2018
- cmake >= 3.9

On Ubuntu (tested with 22.04 LTS) the following should do it:
```bash
sudo apt install cmake g++ zlib1g-dev libtbb-dev libboost-dev libboost-filesystem-dev libboost-program-options-dev
```

### GPU Support (optional)
- NVML or CUDA >= 11.0

Ubuntu (at least 21.10): 
```bash
sudo apt install nvidia-common nvidia-cuda-toolkit nvidia-settings nvidia-driver-XXX nvidia-utils-XXX
```
(Replace XXX with an appropriate driver number, e.g. 470. `nvidia-settings` and `nvidia-utils` are optional but useful.)

### FPGA Support (optional)
- ADMXRC3 (Alpha Data API)
The Alpha Data API is not included in a repository and has to be installed manually.

## Install
Use CMake to configure your project.

### Standard (Release build):
```bash
cmake -DCMAKE_BUILD_TYPE=Release -S src -B build
cd build
make -j
```

### With installed CUDA or NVML but without GPU Support
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DGPU_SUPPORT=OFF -S src -B build
```

### With installed Alpha Data API but without FPGA Support
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DFPGA_SUPPORT=OFF -S src -B build
```
