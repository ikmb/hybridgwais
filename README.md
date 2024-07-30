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
1:4999  1:5000  M0P1    M0P2    159.536 693.048 1.47874e-36     6.5411  0.517871
1:1190  1:3226  N1189   N3225   28.8503 2.96464 7.83682e-08     1.08676 0.202329
1:1830  1:2496  N1829   N2495   27.6975 2.24539 1.42134e-07     0.80888 0.153696
...
```
The meanings of each column are explained in the [Results](#results) section below.

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
1:498   1:499   1:500   M0P3    M0P4    M0P5    0.0826113
1:83    1:301   1:462   N82     N300    N461    0.0309345
1:10    1:22    1:239   N9      N21     N238    0.03081
...
```
The meanings of each column are explained in the [Results](#results) section below.


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

## Results
The generated `scores`-files are always sorted by the score column of the first applied test.
Per default, it contains the first 100,000 results. You can adjust this number with the `-n` switch.

### Result columns

Common columns:
- `POS_A`: genetic position in the format *chr:bp* of variant A
- `POS_B`: genetic position in the format *chr:bp* of variant B
- `POS_C`: genetic position in the format *chr:bp* of variant C (optional)
- `SNPID_A`: identifier string of variant A
- `SNPID_B`: identifier string of variant B
- `SNPID_C`: identifier string of variant C (optional)

Logistic regression:
- `LOGR_CHISQ`: chi-squared score
- `LOGR_OR`: odds-ratio
- `LOGR_P-VAL`: p-value of the chi-squared score assuming a chi-squared distribution with one degree of freedom
- `LOGR_BETA`: beta_3 value after the last iteration
- `LOGR_EPS`: standard error after the last iteration

BOOST:
- `BOOST_CHISQ`: chi-squared score of log-linear test applied after KSA-filtering
- `BOOST_ERR`: standard error after the last iteration of the log-linear test after KSA-filtering
- `BOOST_P-VAL`: p-value of the chi-squared score assuming a chi-squared distribution with four degrees of freedom

Log-linear test:
- `LL_CHISQ`: Chi-squared score
- `LL_ERR`: standard error after the last iteration
- `LL_P-VAL`: p-value of the chi-squared score assuming a chi-squared distribution with four degrees of freedom

Mutual information:
- `MI`: mutual information I(X_1,X_2;Y) or I(X_1,X_2,X_3;Y) respectively

Information gain:
- `IG`: information gain I(X_1;X_2;Y) or I(X_1;X_2;X_3;Y) respectively

Linkage disequilibrium:

Note that the equation to calculate r^2 allows up to three possible solutions, although in most cases they are all equal. For pairwise tests, we report the highest and the lowest solution:
- `R2_H`: largest (highest) solution for r^2
- `R2_L`: smallest (lowest) solution for r^2

For third order tests we report all largest (highest) pairwise r^2 scores:
- `R2_AB_H`: largest (highest) solution for r^2 for variant pair A and B
- `R2_AC_H`: largest (highest) solution for r^2 for variant pair A and C
- `R2_BC_H`: largest (highest) solution for r^2 for variant pair B and C
