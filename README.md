# HybridGWAIS

## Requirements
- Linux operating system (e.g. Ubuntu, this tool was tested with Ubuntu 22.04 LTS)
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

The tool was tested on a system with Tesla P100 and Tesla V100 GPUs.

### FPGA Support (optional)
- ADMXRC3 (Alpha Data API)
The Alpha Data API is not included in a repository and has to be installed manually.

The tool was tested with an Alpha Data ADM-PCIE-8K5 FPGA accelerator.

## Install
Use CMake to configure your project.

### Standard (Release build):
```bash
cmake -DCMAKE_BUILD_TYPE=Release -S src -B build
cd build
make -j
```
The process should not take longer than a few minutes on a standard computer.

### With installed CUDA or NVML but without GPU Support
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DGPU_SUPPORT=OFF -S src -B build
```

### With installed Alpha Data API but without FPGA Support
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DFPGA_SUPPORT=OFF -S src -B build
```

## Example
Small examples can be found in the [example](example) directory. After you successfully compiled the software, navigate to the example directory and launch one of the example scripts.

### 2nd-order (pairwise) interaction test
For a pairwise interaction test example launch:
```bash
cd example
./run2nd
```
Or if you have compiled with GPU support, you can run:
```bash
./run2nd-gpu
```
Both examples launch a second order (pairwise) logistic regression test on an example dataset located in [example/2ndorder](example/2ndorder).
The runtime is only a few seconds on a standard computer, but keep in my mind, that larger input datasets may increase the runtime dramatically as all pairwise combinations of input variants are tested!

The examples create a `scores`-file located in this folder. The top results should look like this:
```
POS A   POS B   SNPID A SNPID B LOGR_CHISQ      LOGR_OR LOGR_P-VAL      LOGR_BETA       LOGR_EPS
0:0     0:0     M0P1    M0P2    159.536 693.048 1.47874e-36     6.5411  0.517871
0:0     0:0     N1189   N3225   28.8503 2.96464 7.83682e-08     1.08676 0.202329
0:0     0:0     N1829   N2495   27.6975 2.24539 1.42134e-07     0.80888 0.153696
...
```
If you want to run the tool directly without a script, you could simply call:
```bash
../build/hybridgwais 2ndorder/2ndorder -m logreg -n 100 --gpu 0
```

### 3rd-order interaction test
Likewise, an example with a third-order test using the information gain measure can be launched:
```bash
cd example
./run3rd
```
Or if you have compiled with GPU support, you can run:
```bash
./run3rd-gpu
```
The dataset for these examples is located in [example/3rdorder](example/3rdorder).
The runtime is only a few seconds on a standard computer, but keep in my mind, that larger input datasets may increase the runtime dramatically as all triple combinations of input variants are tested!

The examples again create a `scores`-file located in this folder. The top results should look like this:
```
POS A   POS B   POS C   SNPID A SNPID B SNPID C IG
0:0     0:0     0:0     M0P3    M0P4    M0P5    0.0826113
0:0     0:0     0:0     N82     N300    N461    0.0309345
0:0     0:0     0:0     N9      N21     N238    0.03081
```

## How to run *HybridGWAIS* with your data
The simplest way to launch an analysis is to call the tool with your `.bed/.bim/.fam` input (provide only the base name without the file ending) and the statistic method(s) you want to use with the `-m` switch.
```bash
hybridgwais my_data -m logreg
```

Please call the tool with the `--help` switch for additional information and user options.
```bash
hybridgwais --help
```

Note: Please add `hybridgwais` to your `PATH` variable if you want to call it without the complete path location. For example, you could add something like this in your `~/.bashrc`:
```bash
export PATH=$PATH:/path/to/hybridgwais/build
```
