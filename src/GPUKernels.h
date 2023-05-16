/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of HybridGWAIS.
 *
 *    HybridGWAIS is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    HybridGWAIS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HybridGWAIS. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GPUKERNELS_H_
#define GPUKERNELS_H_

#ifdef USE_CUDA_GPU
#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>

#include "GPUEngine.h"

using namespace std;

#define checkCUDAError(err) { checkCUDAError_unwrapped((err), __FILE__, __LINE__); }
static void __attribute__((unused)) checkCUDAError_unwrapped(cudaError_t code, const char *file, int line) {
    if(code != cudaSuccess) {
        stringstream ss;
        ss << file << "(" << line << ")";
        string file_and_line;
        ss >> file_and_line;
        throw thrust::system_error(code, thrust::cuda_category(), file_and_line);
    }
}
#else
#include "CUDAHostWrapper.h"
using namespace std;
#endif

#ifdef USE_CUDA_GPU
// host wrapper for copying the constant symbols
__host__ void copyConstantsToDevice(
        const vector<Method> &methods,
        unsigned long numcases,
        unsigned long numsamples,
        unsigned long snpsize,
        unsigned long tablesize,
        unsigned long dbOffset,
        bool gpuonly);
#endif

// SET ID FUNCTIONS

__device__ __host__ void setID2Way(const uint32_t *ctable, DeviceResult<> result);
__device__ __host__ void setID3Way(const uint32_t *ctable, DeviceResult<> result);

// LD (r^2) CALCULATION FUNCTIONS

__device__ __host__ int cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions);

/*
 * Expect allfrequencies as double[7], frequencies are inserted in the following order: f11,f12,f21,f22,f1x,fx1,k,
 * whereby k is n11/n.
 * All frequencies are provided without the het-het share which might have up to three real solutions. The solutions
 * between 0 and k are provided in the array hethet which is expected as double[3].
 * Furthermore, the number of those solutions is returned.
 * The contents of the provided arrays are overwritten.
 */
__device__ __host__ int calcAlleleFrequenciesAndHetHet_2way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *allfrequencies, double *hethet);

__device__ __host__ void calcLD_2way_FromFreq(const double *all_freq, const double *solutions, int solcnt, double *highLD, double *lowLD);
__device__ __host__ void calcLD_2way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *highLD, double *lowLD);
__device__ __host__ void calcPairwiseLD_3way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *ldAB, double *ldAC, double *ldBC);

// APPROXIMATE P-VALUES

// calculate approximate p-value for chi-squared distribution with 1 df
__device__ __host__ inline double approx_p_1df(double chisq);

// calculate approximate p-value for chi-squared distribution with 4 df
__device__ __host__ inline double approx_p_4df(double chisq);


// CONTINGENCY TABLE CREATION

__device__ __host__ void generateContingencyTable2Way(uint32_t *dest, uint32_t *numcases, uint32_t *numctrls, const unsigned char* dbA, const unsigned char* dbB, unsigned snpsize, unsigned casecount, unsigned a, unsigned b);
__device__ __host__ void generateContingencyTable3Way(uint32_t *dest, uint32_t *numcases, uint32_t *numctrls, const unsigned char* dbA, const unsigned char* dbB, const unsigned char* dbC, unsigned snpsize, unsigned casecount, unsigned a, unsigned b, unsigned c);

// CONTINGENCY TABLE DECODING

// decode 3way contingency table
__device__ void decodeContingencyTable3Way(
        const uint16_t *values,
        uint32_t *ctable,
        uint32_t *numcases, uint32_t *numctrls);

// decode 2way contingency table
__device__ void decodeContingencyTable2Way(
        const uint16_t *values,
        uint32_t *ctable,
        uint32_t *numcases, uint32_t *numctrls);

// GENERAL KERNELS
__global__ void Kernel2Way(uint16_t *values, size_t tablesExpected, ResultView<> resultView, int decompResultsFlags);
__global__ void Kernel3Way(uint16_t *values, size_t tablesExpected, ResultView<> resultView, int decompResultsFlags);

// 3WAY MUTUAL INFORMATION KERNEL
__device__ __host__ void MutualInformationKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 2WAY MUTUAL INFORMATION KERNEL
__device__ __host__ void MutualInformationKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 3WAY INFORMATION GAIN KERNEL
__device__ __host__ void InformationGainKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 2WAY INFORMATION GAIN KERNEL
__device__ __host__ void InformationGainKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// KSA filter for 2way BOOST kernel
__device__ __host__ GPUEngine::score_type filterKSA(const uint32_t* cases, const uint32_t* controls, uint32_t numCases, uint32_t numCtrls,
        GPUEngine::score_type invCases, GPUEngine::score_type invCtrls);

// log-linear filter for 2way BOOST and 2way log-linear kernel
__device__ __host__ void filterLog(const uint32_t* cases, const uint32_t* controls, GPUEngine::score_type *result_score, GPUEngine::score_type *result_error, int *result_iterations);

// 2WAY BOOST KERNEL
__device__ __host__ void BOOSTKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 2WAY LOG-LINEAR KERNEL
__device__ __host__ void LogLinearKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// helper functions for logistic regression kernels
__device__ __host__ void fill_zero(double *array, int cnt);
__device__ __host__ void fill_zero(float *array, int cnt);
__device__ __host__ void compute_p_and_v_2way(const GPUEngine::score_type* beta, const uint32_t* nijk, GPUEngine::score_type* pijk, GPUEngine::score_type* vij);
__device__ __host__ void compute_p_and_v_3way(const GPUEngine::score_type* beta, const uint32_t* nijkl, GPUEngine::score_type* pijkl, GPUEngine::score_type* vijk);
__device__ __host__ void compute_gradient_2way(const uint32_t* nijk, const GPUEngine::score_type* pijk, GPUEngine::score_type* grad);
__device__ __host__ void compute_gradient_3way(const uint32_t* nijkl, const GPUEngine::score_type* pijkl, GPUEngine::score_type* grad);
__device__ __host__ void compute_hessian_2way(const GPUEngine::score_type* vij, GPUEngine::score_type* hh);
__device__ __host__ void compute_hessian_3way(const GPUEngine::score_type* vijk, GPUEngine::score_type* hh);
__device__ __host__ void cholesky_decomposition(const GPUEngine::score_type* aa, GPUEngine::score_type* ll, uint32_t param_cnt);
__device__ __host__ void cholesky_decomposition_float(const float* aa, float* ll, uint32_t param_cnt);
__device__ __host__ void solve_linear_system(const GPUEngine::score_type* ll, const GPUEngine::score_type* yy, GPUEngine::score_type* xx, uint32_t param_cnt);
__device__ __host__ void solve_linear_system_float(const float* ll, const float* yy, float* xx, uint32_t param_cnt);

// 2WAY (PLINK) LOGISTIC REGRESSION KERNEL
__device__ __host__ void LogisticRegressionKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 3WAY LOGISTIC REGRESSION KERNEL
__device__ __host__ void LogisticRegressionKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 2WAY LD KERNEL
__device__ __host__ void LDKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

// 3WAY LD KERNEL
__device__ __host__ void LDKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp);

#endif /* GPUKERNELS_H_ */
