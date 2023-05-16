/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of HybridGWAIS.
 *
 *    Includes parts of PLINK 1.90, copyright (C) 2005-2018 Shaun Purcell,
 *    Christopher Chang.
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

#ifdef USE_CUDA_GPU

extern "C" {
#include <cuda_runtime.h>
}

#else

#include "CUDAHostWrapper.h"

#endif

#include "GPUEngine.h"
#include "GPUKernels.h"

// constants
#define PI 3.1415926535897932
// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
// less tolerant version (2^{-35}) for some exact calculations
#define SMALLISH_EPSILON 0.00000000002910383045673370361328125

// ctable decoding
#define TO3(i,j,k) ((i)*9+(j)*3+(k))
#define TO2(i,j)         ((i)*3+(j))

// thresholds
const GPUEngine::score_type cudaKSAThres = 15; // to check if val >= thres, same as 2*val >= 30 as written in the BOOST paper
const int maxLoglinIter = 33; // subject to further testing

// global data information
__constant__ __device__ Method::Type devMethods[Method::Type::TYPE_MAX]; // reserve an array for the order of methods to be processed
__constant__ __device__ unsigned long devNumMethods = 0;
__constant__ __device__ unsigned long devCaseCount = 0;
__constant__ __device__ unsigned long devSampleCount = 0;
__constant__ __device__ unsigned long devSNPSize = 0;
__constant__ __device__ unsigned long devTableSize = 0;
__constant__ __device__ unsigned long devDBOffset = 0;
// flags
__constant__ __device__ int devGPUonly = 0;

#ifdef USE_CUDA_GPU
// host wrapper for copying the constant symbols
__host__ void copyConstantsToDevice(
        const vector<Method> &methods,
        unsigned long numcases,
        unsigned long numsamples,
        unsigned long snpsize,
        unsigned long tablesize,
        unsigned long dbOffset,
        bool gpuonly) {

        Method::Type mts[Method::Type::TYPE_MAX];
        for (int mi = 0; mi < methods.size(); mi++) {
            mts[mi] = methods[mi].getType();
        }
        checkCUDAError(cudaMemcpyToSymbol(devMethods, &mts, sizeof(devMethods), 0, cudaMemcpyHostToDevice))
        unsigned long numMethods = methods.size();
        checkCUDAError(cudaMemcpyToSymbol(devNumMethods, &numMethods, sizeof(devNumMethods), 0, cudaMemcpyHostToDevice))
        checkCUDAError(cudaMemcpyToSymbol(devCaseCount, &numcases, sizeof(devCaseCount), 0, cudaMemcpyHostToDevice))
        checkCUDAError(cudaMemcpyToSymbol(devSampleCount, &numsamples, sizeof(devSampleCount), 0, cudaMemcpyHostToDevice))
        checkCUDAError(cudaMemcpyToSymbol(devSNPSize, &snpsize, sizeof(devSNPSize), 0, cudaMemcpyHostToDevice))
        checkCUDAError(cudaMemcpyToSymbol(devTableSize, &tablesize, sizeof(tablesize), 0, cudaMemcpyHostToDevice))
        checkCUDAError(cudaMemcpyToSymbol(devDBOffset, &dbOffset, sizeof(dbOffset), 0, cudaMemcpyHostToDevice))

        int gpuonlyInt = gpuonly ? 1 : 0;
        checkCUDAError(cudaMemcpyToSymbol(devGPUonly, &gpuonlyInt, sizeof(devGPUonly), 0, cudaMemcpyHostToDevice))
}
#endif


// SET ID FUNCTIONS

__device__ __host__ void setID2Way(const uint32_t *ctable, DeviceResult<> result) {
    const int SNP_A = 1;
    const int SNP_B = 0;

    const uint32_t *id = ctable + 18;

    result.setID(0, id[SNP_A]);
    result.setID(1, id[SNP_B]);
}

__device__ __host__ void setID3Way(const uint32_t *ctable, DeviceResult<> result) {
    const int SNP_A = 2;
    const int SNP_B = 1;
    const int SNP_C = 0;

    const uint32_t *id = ctable + 54;

    result.setID(0, id[SNP_A]);
    result.setID(1, id[SNP_B]);
    result.setID(2, id[SNP_C]);
}

// END SET ID FUNCTIONS


// LD (r^2) CALCULATION FUNCTIONS

__device__ __host__ int cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions) {
  // Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
  // solutions[] (sorted from smallest to largest), and returning the count.
  // Multiple roots are only returned/counted once.
  // Additional research into numerical stability may be in order here.
  double a2 = coef_a * coef_a;
  double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
  double rr = (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
  double r2 = rr * rr;
  double q3 = qq * qq * qq;
  double adiv3 = coef_a * (1.0 / 3.0);
  double sq;
  double dxx;
  if (r2 < q3) {
    // three real roots
    sq = sqrt(qq);
    dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
    sq *= -2;
    solutions[0] = sq * cos(dxx) - adiv3;
    solutions[1] = sq * cos(dxx + (2.0 * PI / 3.0)) - adiv3;
    solutions[2] = sq * cos(dxx - (2.0 * PI / 3.0)) - adiv3;
    // now sort and check for within-epsilon equality
    if (solutions[0] > solutions[1]) {
      dxx = solutions[0];
      solutions[0] = solutions[1];
      if (dxx > solutions[2]) {
        solutions[1] = solutions[2];
    solutions[2] = dxx;
      } else {
    solutions[1] = dxx;
      }
      if (solutions[0] > solutions[1]) {
    dxx = solutions[0];
    solutions[0] = solutions[1];
    solutions[1] = dxx;
      }
    } else if (solutions[1] > solutions[2]) {
      dxx = solutions[1];
      solutions[1] = solutions[2];
      solutions[2] = dxx;
    }
    if (solutions[1] - solutions[0] < EPSILON) {
      solutions[1] = solutions[2];
      return (solutions[1] - solutions[0] < EPSILON)? 1 : 2;
    }
    return (solutions[2] - solutions[1] < EPSILON)? 2 : 3;
  }
  dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
  if (dxx < SMALLISH_EPSILON && dxx > -SMALLISH_EPSILON) { // effectively if dxx == 0
    solutions[0] = -adiv3;
    return 1;
  }
  if (rr < 0.0) {
    dxx = -dxx;
  }
  sq = qq / dxx;
  solutions[0] = dxx + sq - adiv3;
  // use of regular epsilon here has actually burned us
  if (fabs(dxx - sq) >= (EPSILON * 8)) {
    return 1;
  }
  if (dxx >= 0.0) {
    solutions[1] = solutions[0];
    solutions[0] = -dxx - adiv3;
  } else {
    solutions[1] = -dxx - adiv3;
  }
  return 2;
}

/*
 * Expect allfrequencies as double[7], frequencies are inserted in the following order: f11,f12,f21,f22,f1x,fx1,k,
 * whereby k is n11/n.
 * All frequencies are provided without the het-het share which might have up to three real solutions. The solutions
 * between 0 and k are provided in the array hethet which is expected as double[3].
 * Furthermore, the number of those solutions is returned.
 * The contents of the provided arrays are overwritten.
 */
__device__ __host__ int calcAlleleFrequenciesAndHetHet_2way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *allfrequencies, double *hethet) {

    uint64_t f11i = (double)(2*cases[0] + cases[1] + cases[3]);
    uint64_t f12i = (double)(2*cases[2] + cases[1] + cases[5]);
    uint64_t f21i = (double)(2*cases[6] + cases[3] + cases[7]);
    uint64_t f22i = (double)(2*cases[8] + cases[7] + cases[5]);
    f11i += (double)(2*controls[0] + controls[1] + controls[3]);
    f12i += (double)(2*controls[2] + controls[1] + controls[5]);
    f21i += (double)(2*controls[6] + controls[3] + controls[7]);
    f22i += (double)(2*controls[8] + controls[7] + controls[5]);

    uint64_t ki = cases[4] + controls[4];

    uint64_t num_samples = numCases + numControls;
    double divider = 1.0 / (2*num_samples);

    double f11 = f11i * divider;
    double f12 = f12i * divider;
    double f21 = f21i * divider;
    double f22 = f22i * divider;
    double k   = ki * divider;

    double f1x = f11 + f12 + k;
    //double f2x = 1.0 - f1x;
    double fx1 = f11 + f21 + k;
    //double fx2 = 1.0 - fx1;

    int l = 0, r = 0;
    double x[3] = {0.0, 0.0, 0.0};

    // only necessary to divide up K if it is non-zero
    if (ki > 0) {
        // if one of {f11,f22} and one of {f12,f21} is zero this forms a special case
        if (!((f11i == 0 || f22i == 0) && (f12i == 0 || f21i == 0))) { // common case!
            // (f11+x) (f22+x) (K-x) = (f12+(K-x)) (f21+(K-x)) x
            // <=> x^3 + 0.5(f11+f22-f12-f21-3K) x^2 + 0.5(f11f22 + f12f21 + K(f12+f21-f11-f22+K)) x - 0.5(Kf11f22) = 0
            // -> up to three (real) solutions possible
            int solcnt = cubic_real_roots(0.5 * (double)(f11 + f22 - f12 - f21 - 3*k), 0.5 * (double)(f11 * f22 + f12 * f21 + k*(f12 + f21 - f11 - f22 + k)), -0.5 * (double)(k * f11 * f22), x);

            // filter solutions that do not fit (i.e. negative ones and those larger than K)
            // according to PLINK code, at least one solution has to fit,
            // but as we encountered occasional segfaults, we introduce a limit now
            r = solcnt-1;
            // filter negative solutions
            while (l < 2 && x[l] < -SMALLISH_EPSILON)
                l++;
            // set almost zero results to zero
            if (x[l] < SMALLISH_EPSILON)
                x[l] = 0.0;
            // filter solutions > K
            while (r > l && x[r] > k + SMALLISH_EPSILON)
                r--;
            // set results very near to K to exactly K
            if (x[r] > k - SMALLISH_EPSILON)
                x[r] = k;
        } else { // special case with zeros...
            // above equation implies that 0 and K are always solutions, so find only the third one.
            x[0] = 0.0;
            double nz_fxx = f11 + f22;
            double nz_fxy = f12 + f21;
            if ((nz_fxx + SMALLISH_EPSILON < k + nz_fxy) && (nz_fxy + SMALLISH_EPSILON < k + nz_fxx)) {
                r = 2;
                x[1] = (k + nz_fxy - nz_fxx) * 0.5;
                x[2] = k;
            } else {
                r = 1;
                x[1] = k;
            }
        }
    }

    allfrequencies[0] = f11;
    allfrequencies[1] = f12;
    allfrequencies[2] = f21;
    allfrequencies[3] = f22;
    allfrequencies[4] = f1x;
    allfrequencies[5] = fx1;
    allfrequencies[6] = k;

    for (int i=l,j=0; i <= r; i++,j++) {
        hethet[j] = x[i];
    }

    return r-l+1;
}

__device__ __host__ void calcLD_2way_FromFreq(const double *all_freq, const double *solutions, int solcnt, double *highLD, double *lowLD) { //, uint32_t idA, uint32_t idB) {

    const double *x = solutions;

    double f11 = all_freq[0];
//    double f12 = all_freq[1];
//    double f21 = all_freq[2];
//    double f22 = all_freq[3];
    double f1x = all_freq[4];
    double fx1 = all_freq[5];
    double f2x = 1.0 - f1x;
    double fx2 = 1.0 - fx1;
//    double k = all_freq[6];

    // Solutions are ordered and r^2 is a positive quadratic function, i.e.
    // the leftmost solution leads either to the smallest or largest LD,
    // then the rightmost is either the largest or the smallest respectively.
    // In the case of three solutions, the one in the middle will neither be largest nor smallest.

    // solution 1
    // D = f11+x - fx1*f1x
    double d = f11 + x[0] - fx1*f1x;
    // r^2 = D^2 / (f1x*fx1*f2x*fx2)
    *lowLD = d*d / (f1x*fx1*f2x*fx2); // assuming smallest

    // solution 2 (if different to 1)
    if (solcnt > 1) {
        d = f11 + x[solcnt-1] - fx1*f1x;
        double rsq = d*d / (f1x*fx1*f2x*fx2);
        if (rsq < *lowLD) {
            *highLD = *lowLD;
            *lowLD = rsq;
        } else
            *highLD = rsq;
    } else
        *highLD = *lowLD;

}

__device__ __host__ void calcLD_2way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *highLD, double *lowLD) { //, uint32_t idA, uint32_t idB) {
    double allfreq[7], x[3];
    int solcnt = calcAlleleFrequenciesAndHetHet_2way(cases, controls, numCases, numControls, allfreq, x);
    calcLD_2way_FromFreq(allfreq, x, solcnt, highLD, lowLD); //, idA, idB);
}

__device__ __host__ void calcPairwiseLD_3way(const uint32_t *cases, const uint32_t *controls, uint32_t numCases, uint32_t numControls, double *ldAB, double *ldAC, double *ldBC) { //, uint32_t idA, uint32_t idB, uint32_t idC) {
    uint32_t tmpcases[9], tmpctrls[9];
    double rsqh, rsql;

    // pair AB
    for (int i = 0, j = 0; i < 9; i++, j+=3) {
        tmpcases[i] = cases[j] + cases[j+1] + cases[j+2];
        tmpctrls[i] = controls[j] + controls[j+1] + controls[j+2];
    }
    calcLD_2way(tmpcases, tmpctrls, numCases, numControls, &rsqh, &rsql); //, idA, idB);
    *ldAB = rsqh;

    // pair AC
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int k = TO3(i,0,j);
            tmpcases[TO2(i,j)] = cases[k] + cases[k+3] + cases[k+6];
            tmpctrls[TO2(i,j)] = controls[k] + controls[k+3] + controls[k+6];
        }
    }
    calcLD_2way(tmpcases, tmpctrls, numCases, numControls, &rsqh, &rsql); //, idA, idC);
    *ldAC = rsqh;

    // pair BC
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int k = TO3(0,i,j);
            tmpcases[TO2(i,j)] = cases[k] + cases[k+9] + cases[k+18];
            tmpctrls[TO2(i,j)] = controls[k] + controls[k+9] + controls[k+18];
        }
    }
    calcLD_2way(tmpcases, tmpctrls, numCases, numControls, &rsqh, &rsql); //, idB, idC);
    *ldBC = rsqh;
}

// END OF LD CALCULATION FUNCTIONS

// APPROXIMATE P-VALUES

// calculate approximate p-value for chi-squared distribution with 1 df (from PLINK)
__device__ __host__ inline double approx_p_1df(double chisq) {
    // original PLINK implementation
    //double zz = -sqrt(chisq);
    //// copied from PLINK normdist()
    //double sqrt2pi = 2.50662827463;
    //double t0;
    //double z1;
    //double p0;
    //t0 = 1 / (1 + 0.2316419 * fabs(zz));
    //z1 = exp(-0.5 * zz * zz) / sqrt2pi;
    //p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
    //// return zz >= 0 ? 1 - p0 : p0; -> in this case zz is always negative (or zero, but then the result doesn't matter anyway)
    //return p0 * 2; // *2 is what PLINK does after the function call anyway, so I put in here.

    // slightly improved
    double sqrt2pi = 2.50662827463;
    double t0 = 1 / (1 + 0.2316419 * sqrt(chisq));
    double z1 = exp(-0.5 * chisq) / sqrt2pi;
    double p0 = z1 * t0 * (0.31938153 + t0 * (-0.356563782 + t0 * (1.781477937 + t0 * (-1.821255978 + 1.330274429 * t0))));
    return p0 * 2; // *2 is what PLINK does after the function call anyway, so I put in here.
}

// calculate approximate p-value for chi-squared distribution with 4 df
__device__ __host__ inline double approx_p_4df(double chisq) {
    double chisqhalf = 0.5 * chisq;
    return (1 + chisqhalf) * exp(-chisqhalf);
}

//// calculate approximate p-value for chi-squared distribution with 8 df
//__device__ inline double approx_p_8df(double chisq) {
//    double chisqhalf = 0.5 * chisq;
//    return (1 + chisqhalf * (1 + chisqhalf * (0.5 + chisqhalf/6))) * exp(-chisqhalf);
//}

// END OF APPROXIMATE P-VALUES


// CONTINGENCY TABLE CREATION

__device__ __host__ void generateContingencyTable2Way(uint32_t *dest, uint32_t *numcases, uint32_t *numctrls,
        const unsigned char* dbA, const unsigned char* dbB, unsigned snpsize, unsigned casecount, unsigned a, unsigned b) {

    // initialize counters
    for (int i = 0; i < 18; i++)
        dest[i] = 0;

    unsigned gtpos = 0;

    for(unsigned pos = 0; pos < snpsize; pos++) {

        unsigned char gta = dbA[pos] & 3;
        unsigned char gtb = dbB[pos] & 3;
        unsigned counter = gta * 3 + gtb;
        if(gtpos < casecount)
            counter += 9;
        if (gta != 3 && gtb != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;

        gta = (dbA[pos] >> 2) & 3;
        gtb = (dbB[pos] >> 2) & 3;
        counter = gta * 3 + gtb;
        if(gtpos < casecount)
            counter += 9;
        if (gta != 3 && gtb != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;

        gta = (dbA[pos] >> 4) & 3;
        gtb = (dbB[pos] >> 4) & 3;
        counter = gta * 3 + gtb;
        if(gtpos < casecount)
            counter += 9;
        if (gta != 3 && gtb != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;

        gta = (dbA[pos] >> 6) & 3;
        gtb = (dbB[pos] >> 6) & 3;
        counter = gta * 3 + gtb;
        if(gtpos < casecount)
            counter += 9;
        if (gta != 3 && gtb != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;
    };

    // ID
    dest[18] = b;
    dest[19] = a;

    *numctrls = 0;
    for (int i = 0; i < 9; i++)
        *numctrls += dest[i];

    *numcases = 0;
    for (int i = 9; i < 18; i++)
        *numcases += dest[i];
}

__device__ __host__ void generateContingencyTable3Way(uint32_t *dest, uint32_t *numcases, uint32_t *numctrls,
        const unsigned char* dbA, const unsigned char* dbB, const unsigned char* dbC, unsigned snpsize, unsigned casecount, unsigned a, unsigned b, unsigned c) {

    // initialize counters
    for (int i = 0; i < 54; i++)
        dest[i] = 0;

    unsigned gtpos = 0;

    for(unsigned pos = 0; pos < snpsize; pos++) {
        unsigned char gta = dbA[pos] & 3;
        unsigned char gtb = dbB[pos] & 3;
        unsigned char gtc = dbC[pos] & 3;
        unsigned counter = gta * 9 + gtb * 3 + gtc;

        if(gtpos < casecount)
            counter += 27;
        if (gta != 3 && gtb != 3 && gtc != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;

        gta = (dbA[pos] >> 2) & 3;
        gtb = (dbB[pos] >> 2) & 3;
        gtc = (dbC[pos] >> 2) & 3;
        counter = gta * 9 + gtb * 3 + gtc;
        if(gtpos < casecount)
            counter += 27;
        if (gta != 3 && gtb != 3 && gtc != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;
        gta = (dbA[pos] >> 4) & 3;
        gtb = (dbB[pos] >> 4) & 3;
        gtc = (dbC[pos] >> 4) & 3;
        counter = gta * 9 + gtb * 3 + gtc;
        if(gtpos < casecount)
            counter += 27;
        if (gta != 3 && gtb != 3 && gtc != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;
        gta = (dbA[pos] >> 6) & 3;
        gtb = (dbB[pos] >> 6) & 3;
        gtc = (dbC[pos] >> 6) & 3;
        counter = gta * 9 + gtb * 3 + gtc;
        if(gtpos < casecount)
            counter += 27;
        if (gta != 3 && gtb != 3 && gtc != 3) // do not count unknowns
            dest[counter]++;
        gtpos++;
    };

    // ID
    dest[54] = c;
    dest[55] = b;
    dest[56] = a;

    *numctrls = 0;
    for (int i = 0; i < 27; i++)
        *numctrls += dest[i];

    *numcases = 0;
    for (int i = 27; i < 54; i++)
        *numcases += dest[i];

}

// END OF CONTINGENCY TABLE CREATION


// CONTINGENCY TABLE DECODING

// decode 3way contingency table
__device__ void decodeContingencyTable3Way(
        const uint16_t *values,
        uint32_t *ctable,
        uint32_t *numcases, uint32_t *numctrls) {

    const int SNP_A = 2;
    const int SNP_B = 1;
    const int SNP_C = 0;

    uint32_t *controls = ctable;
    uint32_t *cases = &ctable[27];
    uint32_t *id = &ctable[54];

    // decode 16 bit packed values to regular packed 32 bit values, so no decoding required, just mapping

    *numctrls = 0;
    for (int i = 0; i < 27; i++)
        *numctrls += values[i];

    *numcases = 0;
    for (int i = 27; i < 54; i++)
        *numcases += values[i];

    controls[ 0] = values[ 0];
    controls[ 1] = values[ 1];
    controls[ 2] = values[ 2];
    controls[ 3] = values[ 3];
    controls[ 4] = values[ 4];
    controls[ 5] = values[ 5];
    controls[ 6] = values[ 6];
    controls[ 7] = values[ 7];
    controls[ 8] = values[ 8];
    controls[ 9] = values[ 9];
    controls[10] = values[10];
    controls[11] = values[11];
    controls[12] = values[12];
    controls[13] = values[13];
    controls[14] = values[14];
    controls[15] = values[15];
    controls[16] = values[16];
    controls[17] = values[17];
    controls[18] = values[18];
    controls[19] = values[19];
    controls[20] = values[20];
    controls[21] = values[21];
    controls[22] = values[22];
    controls[23] = values[23];
    controls[24] = values[24];
    controls[25] = values[25];
    controls[26] = values[26];
    cases[ 0] = values[27];
    cases[ 1] = values[28];
    cases[ 2] = values[29];
    cases[ 3] = values[30];
    cases[ 4] = values[31];
    cases[ 5] = values[32];
    cases[ 6] = values[33];
    cases[ 7] = values[34];
    cases[ 8] = values[35];
    cases[ 9] = values[36];
    cases[10] = values[37];
    cases[11] = values[38];
    cases[12] = values[39];
    cases[13] = values[40];
    cases[14] = values[41];
    cases[15] = values[42];
    cases[16] = values[43];
    cases[17] = values[44];
    cases[18] = values[45];
    cases[19] = values[46];
    cases[20] = values[47];
    cases[21] = values[48];
    cases[22] = values[49];
    cases[23] = values[50];
    cases[24] = values[51];
    cases[25] = values[52];
    cases[26] = values[53];

    id[SNP_C] =  (static_cast<uint32_t>(values[54]) | static_cast<uint32_t>(values[55]) << 16) + 1; // ID C is off-by-one
    id[SNP_B] =  (static_cast<uint32_t>(values[56]) | static_cast<uint32_t>(values[57]) << 16);     // ID B
    id[SNP_A] =  (static_cast<uint32_t>(values[58]) | static_cast<uint32_t>(values[59]) << 16);     // ID A

}

// decode 2way contingency table
__device__ void decodeContingencyTable2Way(
        const uint16_t *values,
        uint32_t *ctable,
        uint32_t *numcases, uint32_t *numctrls) {

    const int SNP_A = 1;
    const int SNP_B = 0;

    uint32_t *controls = ctable;
    uint32_t *cases = &ctable[9];
    uint32_t *id = &ctable[18];

    // decode 16 bit packed values to regular packed 32 bit values, so no decoding required, just mapping

    *numctrls = 0;
    for (int i = 0; i < 9; i++)
        *numctrls += values[i];

    *numcases = 0;
    for (int i = 9; i < 18; i++)
        *numcases += values[i];

    controls[0] = values[0];
    controls[1] = values[1];
    controls[2] = values[2];
    controls[3] = values[3];
    controls[4] = values[4];
    controls[5] = values[5];
    controls[6] = values[6];
    controls[7] = values[7];
    controls[8] = values[8];
    cases[0] = values[9];
    cases[1] = values[10];
    cases[2] = values[11];
    cases[3] = values[12];
    cases[4] = values[13];
    cases[5] = values[14];
    cases[6] = values[15];
    cases[7] = values[16];
    cases[8] = values[17];

    id[SNP_B] =  (static_cast<uint32_t>(values[18]) | static_cast<uint32_t>(values[19]) << 16) + 1; // ID B is off-by-one
    id[SNP_A] =  (static_cast<uint32_t>(values[20]) | static_cast<uint32_t>(values[21]) << 16);    // ID A

}

// END OF CONTINGENCY TABLE DECODING

// GENERAL KERNEL ACTIVATION FUNCTIONS

__global__ void Kernel2Way(uint16_t *values, size_t tablesExpected, ResultView<> resultView, int decompResultsFlags_) {
    const unsigned char* db = (unsigned char*)(values + devDBOffset/sizeof(*values)); // only required for GPU only run

    uint32_t gid = threadIdx.x + blockIdx.x * blockDim.x; /*global id*/
    values += gid  * (devTableSize/sizeof(*values));

    if(gid >= tablesExpected)
        return;

    uint32_t ctable[20];

    uint32_t numCases = 0;
    uint32_t numControls = 0;

    if (devGPUonly) {
        unsigned a = ((unsigned*)values)[0];
        unsigned b = ((unsigned*)values)[1];
        generateContingencyTable2Way(ctable, &numCases, &numControls,
                db + a*devSNPSize,
                db + b*devSNPSize,
                devSNPSize, devCaseCount, a, b);
    } else
        decodeContingencyTable2Way(values, ctable, &numCases, &numControls);

    setID2Way(ctable, resultView.getDeviceResult(gid));

    int nextResultIndex = 0;
    int decompResultsFlags = decompResultsFlags_;
    for (size_t mi = 0; mi < devNumMethods; mi++) {
        switch (devMethods[mi]) {
        case Method::Type::MI:
            MutualInformationKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::IG:
            InformationGainKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::BOOST:
            BOOSTKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::LOGLIN:
            LogLinearKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::LOGREG:
            LogisticRegressionKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::LD:
            LDKernel2WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        default: // 3-way or unknown method should have been catched already
            break;
        }
        decompResultsFlags >>= 1; // next flag
    }
}

__global__ void Kernel3Way(uint16_t *values, size_t tablesExpected, ResultView<> resultView, int decompResultsFlags_) {
    const unsigned char* db = (unsigned char*)(values + devDBOffset/sizeof(*values)); // only required for GPU only run

    uint32_t gid = threadIdx.x + blockIdx.x * blockDim.x; /*global id*/
    values += gid  * (devTableSize/sizeof(*values));

    if(gid >= tablesExpected)
        return;

    uint32_t ctable[57];

    uint32_t numCases = 0;
    uint32_t numControls = 0;

    if (devGPUonly) {
        unsigned a = ((unsigned*)values)[0];
        unsigned b = ((unsigned*)values)[1];
        unsigned c = ((unsigned*)values)[2];
        generateContingencyTable3Way(ctable, &numCases, &numControls,
                db + a*devSNPSize,
                db + b*devSNPSize,
                db + c*devSNPSize,
                devSNPSize, devCaseCount, a, b, c);
    } else
        decodeContingencyTable3Way(values, ctable, &numCases, &numControls);

    setID3Way(ctable, resultView.getDeviceResult(gid));

    int nextResultIndex = 0;
    int decompResultsFlags = decompResultsFlags_;
    for (size_t mi = 0; mi < devNumMethods; mi++) {
        switch (devMethods[mi]) {
        case Method::Type::MI3:
            MutualInformationKernel3WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::IG3:
            InformationGainKernel3WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::LOGREG3:
            LogisticRegressionKernel3WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        case Method::Type::LD3:
            LDKernel3WayCore(ctable, numCases, numControls, resultView.getDeviceResult(gid), nextResultIndex, decompResultsFlags & 0x1);
            break;
        default: // 2-way or unknown method should have been catched already
            break;
        }
        decompResultsFlags >>= 1; // next flag
    }
}

// END GENERAL KERNEL ACTIVATION FUNCTIONS



// 3WAY MUTUAL INFORMATION KERNEL

__device__ __host__ void MutualInformationKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 27;

    uint32_t no_samples = numCases + numControls;

    /////////////////////////////////////////////////////////////
    // calculate mutual information I(X1,X2,X3;Y)
    /////////////////////////////////////////////////////////////

    GPUEngine::score_type score(0);
    GPUEngine::score_type hxycase(0);
    GPUEngine::score_type hxyctrl(0);
    GPUEngine::score_type hx(0);

    GPUEngine::score_type hy_pre = - (numCases    * log(static_cast<GPUEngine::score_type>(numCases)))
                                     - (numControls * log(static_cast<GPUEngine::score_type>(numControls)));

    uint32_t case_val;
    uint32_t ctrl_val;

    for(int i=0; i < 3; ++i)
        for(int j=0; j < 3; ++j)
            for(int k=0; k < 3; ++k) {
                case_val = cases[TO3(i,j,k)];
                ctrl_val = controls[TO3(i,j,k)];

                hxycase += (case_val == 0 ? 0 : case_val * log(static_cast<GPUEngine::score_type>(case_val)));
                hxyctrl += (ctrl_val == 0 ? 0 : ctrl_val * log(static_cast<GPUEngine::score_type>(ctrl_val)));
                hx += ((case_val + ctrl_val) == 0 ? 0 : (case_val + ctrl_val) * log(static_cast<GPUEngine::score_type>(case_val + ctrl_val)));
            }

    score = hxycase + hxyctrl - hx;

    // backtransform score
    score += hy_pre;
    score /= no_samples;
    score += log(static_cast<GPUEngine::score_type>(no_samples));

    // write back results
    result.setScore(nextResultIndex++, score);

    if (decomp) {
        GPUEngine::score_type hxy = -hxycase - hxyctrl;
        hx = -hx;
        hxy /= no_samples;
        hx /= no_samples;
        hxy += log(static_cast<GPUEngine::score_type>(no_samples));
        hx += log(static_cast<GPUEngine::score_type>(no_samples));
        result.setScore(nextResultIndex++, hxy);
        result.setScore(nextResultIndex++, hx);
    }
}

// END OF 3WAY MUTUAL INFORMATION KERNEL


// 2WAY MUTUAL INFORMATION KERNEL

__device__ __host__ void MutualInformationKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 9;

    uint32_t no_samples = numCases + numControls;

    /////////////////////////////////////////////////////////////
    // calculate mutual information I(X1,X2;Y)
    /////////////////////////////////////////////////////////////

    GPUEngine::score_type score(0);
    GPUEngine::score_type hxycase(0);
    GPUEngine::score_type hxyctrl(0);
    GPUEngine::score_type hx(0);

    GPUEngine::score_type hy_pre = - (numCases    * log(static_cast<GPUEngine::score_type>(numCases)))
                                     - (numControls * log(static_cast<GPUEngine::score_type>(numControls)));

    uint32_t case_val;
    uint32_t ctrl_val;

    for(int i=0; i < 3; ++i)
        for(int j=0; j < 3; ++j) {
            case_val = cases[TO2(i,j)];
            ctrl_val = controls[TO2(i,j)];

            hxycase     += (case_val == 0 ? 0 : case_val * log(static_cast<GPUEngine::score_type>(case_val)));
            hxyctrl     += (ctrl_val == 0 ? 0 : ctrl_val * log(static_cast<GPUEngine::score_type>(ctrl_val)));
            hx += ((case_val + ctrl_val) == 0 ? 0 : (case_val + ctrl_val) * log(static_cast<GPUEngine::score_type>(case_val + ctrl_val)));

        }

    score = hxycase + hxyctrl - hx;

    // backtransform score
    score += hy_pre;
    score /= no_samples;
    score += log(static_cast<GPUEngine::score_type>(no_samples));

    // write back results
    result.setScore(nextResultIndex++, score);

    if (decomp) {
        GPUEngine::score_type hxy =  -hxycase - hxyctrl;
        hx = -hx;
        hxy /= no_samples;
        hx /= no_samples;
        hxy += log(static_cast<GPUEngine::score_type>(no_samples));
        hx += log(static_cast<GPUEngine::score_type>(no_samples));
        result.setScore(nextResultIndex++, hxy);
        result.setScore(nextResultIndex++, hx);
    }
}

// END OF 2WAY MUTUAL INFORMATION KERNEL



// 3WAY INFORMATION GAIN KERNEL

__device__ __host__ void InformationGainKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 27;

    uint32_t no_samples = numCases + numControls;

    uint32_t ac_snpA_controls[3] = {0,0,0};
    uint32_t ac_snpA_cases[3] = {0,0,0};
    uint32_t ac_snpB_controls[3] = {0,0,0};
    uint32_t ac_snpB_cases[3] = {0,0,0};
    uint32_t ac_snpC_controls[3] = {0,0,0};
    uint32_t ac_snpC_cases[3] = {0,0,0};
    for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				ac_snpA_controls[i] += controls[TO3(i,j,k)];
				ac_snpA_cases[i] += cases[TO3(i,j,k)];
				ac_snpB_controls[j] += controls[TO3(i,j,k)];
				ac_snpB_cases[j] += cases[TO3(i,j,k)];
				ac_snpC_controls[k] += controls[TO3(i,j,k)];
				ac_snpC_cases[k] += cases[TO3(i,j,k)];
			}
		}
	}

    /////////////////////////////////////////////////////////////
    // calculate information gain I(X1;X2;X3;Y)
    /////////////////////////////////////////////////////////////

    GPUEngine::score_type   score(0);
    GPUEngine::score_type   mi3(0);
    GPUEngine::score_type   mi2a(0);
    GPUEngine::score_type   mi2b(0);
    GPUEngine::score_type   mi2c(0);
    GPUEngine::score_type   mi1a(0);
    GPUEngine::score_type   mi1b(0);
    GPUEngine::score_type   mi1c(0);

    GPUEngine::score_type hy_pre = - (numCases    * log(static_cast<GPUEngine::score_type>(numCases)))
                                     - (numControls * log(static_cast<GPUEngine::score_type>(numControls)));

    uint32_t case_val;
    uint32_t ctrl_val;
    uint32_t sum_i_cases[9] = {0,0,0,0,0,0,0,0,0};
    uint32_t sum_i_ctrls[9] = {0,0,0,0,0,0,0,0,0};
    uint32_t sum_j_cases[3] = {0,0,0};
    uint32_t sum_j_ctrls[3] = {0,0,0};
    uint32_t sum_k_cases = 0;
    uint32_t sum_k_ctrls = 0;
    GPUEngine::score_type tmp;

    for(int i=0; i < 3; ++i) {
        // the single entropies H(X1), H(X2), H(X3)
        tmp = (ac_snpA_controls[i]+ac_snpA_cases[i]) * log(static_cast<GPUEngine::score_type>(ac_snpA_controls[i]+ac_snpA_cases[i]));
        mi1a -= (::isnan(tmp) ? 0 : tmp);
        tmp = (ac_snpB_controls[i]+ac_snpB_cases[i]) * log(static_cast<GPUEngine::score_type>(ac_snpB_controls[i]+ac_snpB_cases[i]));
        mi1b -= (::isnan(tmp) ? 0 : tmp);
        tmp = (ac_snpC_controls[i]+ac_snpC_cases[i]) * log(static_cast<GPUEngine::score_type>(ac_snpC_controls[i]+ac_snpC_cases[i]));
        mi1c -= (::isnan(tmp) ? 0 : tmp);

        // the pairwise entropies H(X1,Y), H(X2,Y), H(X3,Y)
        tmp = ac_snpA_controls[i] * log(static_cast<GPUEngine::score_type>(ac_snpA_controls[i]));
        mi1a += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpB_controls[i] * log(static_cast<GPUEngine::score_type>(ac_snpB_controls[i]));
        mi1b += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpC_controls[i] * log(static_cast<GPUEngine::score_type>(ac_snpC_controls[i]));
        mi1c += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpA_cases[i] * log(static_cast<GPUEngine::score_type>(ac_snpA_cases[i]));
        mi1a += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpB_cases[i] * log(static_cast<GPUEngine::score_type>(ac_snpB_cases[i]));
        mi1b += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpC_cases[i] * log(static_cast<GPUEngine::score_type>(ac_snpC_cases[i]));
        mi1c += (::isnan(tmp) ? 0 : tmp);

        for (int k=0; k < 3; ++k) {
            sum_j_cases[k] = 0;
            sum_j_ctrls[k] = 0;
        }

        for(int j=0; j < 3; ++j) {

            sum_k_cases = 0;
            sum_k_ctrls = 0;

            for(int k=0; k < 3; ++k) {

                case_val = cases[TO3(i,j,k)];
                ctrl_val = controls[TO3(i,j,k)];

                sum_i_cases[j*3+k] += case_val;
                sum_i_ctrls[j*3+k] += ctrl_val;
                sum_j_cases[k] += case_val;
                sum_j_ctrls[k] += ctrl_val;
                sum_k_cases += case_val;
                sum_k_ctrls += ctrl_val;

                // the combined entropies H(X1,X2,X3,Y) and H(X1,X2,X3)
                tmp = case_val * log(static_cast<GPUEngine::score_type>(case_val));
                mi3 += (::isnan(tmp) ? 0 : tmp);
                tmp = ctrl_val * log(static_cast<GPUEngine::score_type>(ctrl_val));
                mi3 += (::isnan(tmp) ? 0 : tmp);
                tmp = (case_val + ctrl_val) * log(static_cast<GPUEngine::score_type>(case_val + ctrl_val));
                mi3 -= (::isnan(tmp) ? 0 : tmp);

            } //k

            // entropies H(X1,X2,Y) and H(X1,X2)
            tmp = sum_k_cases * log(static_cast<GPUEngine::score_type>(sum_k_cases));
            mi2a += (::isnan(tmp) ? 0 : tmp);
            tmp = sum_k_ctrls * log(static_cast<GPUEngine::score_type>(sum_k_ctrls));
            mi2a += (::isnan(tmp) ? 0 : tmp);
            tmp = (sum_k_cases + sum_k_ctrls) * log(static_cast<GPUEngine::score_type>(sum_k_cases + sum_k_ctrls));
            mi2a -= (::isnan(tmp) ? 0 : tmp);

        } //j

        // entropies H(X1,X3,Y) and H(X1,X3)
        for(int k=0; k < 3; ++k) {
            tmp = sum_j_cases[k] * log(static_cast<GPUEngine::score_type>(sum_j_cases[k]));
            mi2b += (::isnan(tmp) ? 0 : tmp);
            tmp = sum_j_ctrls[k] * log(static_cast<GPUEngine::score_type>(sum_j_ctrls[k]));
            mi2b += (::isnan(tmp) ? 0 : tmp);
            tmp = (sum_j_cases[k] + sum_j_ctrls[k]) * log(static_cast<GPUEngine::score_type>(sum_j_cases[k] + sum_j_ctrls[k]));
            mi2b -= (::isnan(tmp) ? 0 : tmp);
        }

    } //i

    // entropies H(X2,X3,Y) and H(X2,X3)
    for(int jk=0; jk < 9; ++jk) {
        tmp = sum_i_cases[jk] * log(static_cast<GPUEngine::score_type>(sum_i_cases[jk]));
        mi2c += (::isnan(tmp) ? 0 : tmp);
        tmp = sum_i_ctrls[jk] * log(static_cast<GPUEngine::score_type>(sum_i_ctrls[jk]));
        mi2c += (::isnan(tmp) ? 0 : tmp);
        tmp = (sum_i_cases[jk] + sum_i_ctrls[jk]) * log(static_cast<GPUEngine::score_type>(sum_i_cases[jk] + sum_i_ctrls[jk]));
        mi2c -= (::isnan(tmp) ? 0 : tmp);
    }

    score = mi3 - mi2a - mi2b - mi2c + mi1a + mi1b + mi1c;

    // backtransform scores
    score += hy_pre;
    score /= no_samples;
    score += log(static_cast<GPUEngine::score_type>(no_samples));

    // write back results
    result.setScore(nextResultIndex++, score);

    if (decomp) {
//        assert(result.getScoreFields() >= 8);

        // backtransform partial scores
        mi3 += hy_pre;
        mi3 /= no_samples;
        mi3 += log(static_cast<GPUEngine::score_type>(no_samples));
        mi2a += hy_pre;
        mi2a /= no_samples;
        mi2a += log(static_cast<GPUEngine::score_type>(no_samples));
        mi2b += hy_pre;
        mi2b /= no_samples;
        mi2b += log(static_cast<GPUEngine::score_type>(no_samples));
        mi2c += hy_pre;
        mi2c /= no_samples;
        mi2c += log(static_cast<GPUEngine::score_type>(no_samples));
        mi1a += hy_pre;
        mi1a /= no_samples;
        mi1a += log(static_cast<GPUEngine::score_type>(no_samples));
        mi1b += hy_pre;
        mi1b /= no_samples;
        mi1b += log(static_cast<GPUEngine::score_type>(no_samples));
        mi1c += hy_pre;
        mi1c /= no_samples;
        mi1c += log(static_cast<GPUEngine::score_type>(no_samples));

        result.setScore(nextResultIndex++, mi3);
        result.setScore(nextResultIndex++, mi2a);
        result.setScore(nextResultIndex++, mi2b);
        result.setScore(nextResultIndex++, mi2c);
        result.setScore(nextResultIndex++, mi1a);
        result.setScore(nextResultIndex++, mi1b);
        result.setScore(nextResultIndex++, mi1c);
    }
}

// END OF 3WAY INFORMATION GAIN KERNEL


// 2WAY INFORMATION GAIN KERNEL

__device__ __host__ void InformationGainKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 9;

    uint32_t no_samples = numCases + numControls;

    uint32_t ac_snpA_controls[3] = {0,0,0};
    uint32_t ac_snpA_cases[3] = {0,0,0};
    uint32_t ac_snpB_controls[3] = {0,0,0};
    uint32_t ac_snpB_cases[3] = {0,0,0};
    for (int i = 0; i < 3; ++i) {
    	for (int j = 0; j < 3; j++) {
    		ac_snpA_controls[i] += controls[TO2(i,j)];
    		ac_snpA_cases[i] += cases[TO2(i,j)];
    		ac_snpB_controls[j] += controls[TO2(i,j)];
			ac_snpB_cases[j] += cases[TO2(i,j)];
    	}
    }

    /////////////////////////////////////////////////////////////
    // calculate information gain I(X1;X2;Y)
    /////////////////////////////////////////////////////////////

    GPUEngine::score_type   score(0);
    GPUEngine::score_type   mi2(0);
    GPUEngine::score_type   mi1a(0);
    GPUEngine::score_type   mi1b(0);

    GPUEngine::score_type hy_pre = - (numCases    * log(static_cast<GPUEngine::score_type>(numCases)))
                                     - (numControls * log(static_cast<GPUEngine::score_type>(numControls)));

    uint32_t case_val;
    uint32_t ctrl_val;
    GPUEngine::score_type tmp;

    for(int i=0; i < 3; ++i) {
        // the single entropies H(X1), H(X2)
        tmp = (ac_snpA_controls[i]+ac_snpA_cases[i]) * log(static_cast<GPUEngine::score_type>(ac_snpA_controls[i]+ac_snpA_cases[i]));
        mi1a -= (::isnan(tmp) ? 0 : tmp);
        tmp = (ac_snpB_controls[i]+ac_snpB_cases[i]) * log(static_cast<GPUEngine::score_type>(ac_snpB_controls[i]+ac_snpB_cases[i]));
        mi1b -= (::isnan(tmp) ? 0 : tmp);

        // the pairwise entropies H(X1,Y), H(X2,Y)
        tmp = ac_snpA_controls[i] * log(static_cast<GPUEngine::score_type>(ac_snpA_controls[i]));
        mi1a += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpB_controls[i] * log(static_cast<GPUEngine::score_type>(ac_snpB_controls[i]));
        mi1b += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpA_cases[i] * log(static_cast<GPUEngine::score_type>(ac_snpA_cases[i]));
        mi1a += (::isnan(tmp) ? 0 : tmp);
        tmp = ac_snpB_cases[i] * log(static_cast<GPUEngine::score_type>(ac_snpB_cases[i]));
        mi1b += (::isnan(tmp) ? 0 : tmp);

        for(int j=0; j < 3; ++j) {

            case_val = cases[TO2(i,j)];
            ctrl_val = controls[TO2(i,j)];

            // the combined entropies H(X1,X2,Y) and H(X1,X2)
            tmp = case_val * log(static_cast<GPUEngine::score_type>(case_val));
            mi2 += (::isnan(tmp) ? 0 : tmp);
            tmp = ctrl_val * log(static_cast<GPUEngine::score_type>(ctrl_val));
            mi2 += (::isnan(tmp) ? 0 : tmp);
            tmp = (case_val + ctrl_val) * log(static_cast<GPUEngine::score_type>(case_val + ctrl_val));
            mi2 -= (::isnan(tmp) ? 0 : tmp);

        } //j

    } //i

    score = mi2 - mi1a - mi1b;

    // backtransform scores
    // (here we need to remove H(Y) instead of adding it)
    score -= hy_pre;
    score /= no_samples;
    score -= log(static_cast<GPUEngine::score_type>(no_samples));

    // write back results
    result.setScore(nextResultIndex++, score);

    if (decomp) {
//        assert(result.getScoreFields() >= 4);

        // backtransform partial scores
        mi2  += hy_pre;
        mi2  /= no_samples;
        mi2  += log(static_cast<GPUEngine::score_type>(no_samples));
        mi1a += hy_pre;
        mi1a /= no_samples;
        mi1a += log(static_cast<GPUEngine::score_type>(no_samples));
        mi1b += hy_pre;
        mi1b /= no_samples;
        mi1b += log(static_cast<GPUEngine::score_type>(no_samples));

        result.setScore(nextResultIndex++, mi2);
        result.setScore(nextResultIndex++, mi1a);
        result.setScore(nextResultIndex++, mi1b);
    }
}

// END OF 2WAY INFORMATION GAIN KERNEL



// KSA filter for 2way BOOST kernel

__device__ __host__ GPUEngine::score_type filterKSA(const uint32_t* cases, const uint32_t* controls, uint32_t numCases, uint32_t numCtrls,
        GPUEngine::score_type invCases, GPUEngine::score_type invCtrls){
    uint32_t numInds = numCases+numCtrls;

    // Additional for Kirkwood
    // Store the total number for each possible SNP1 and SNP2
    // The order is pi_x0x, pi_x1x, pi_x2x, pi_xx0, pi_xx1, pi_xx2
    uint32_t pi_jk[6];
    // Store the values of p^K_ijk
    GPUEngine::score_type p[18];

    uint32_t ac_snpA_controls[3] = {0,0,0};
    uint32_t ac_snpA_cases[3] = {0,0,0};
    uint32_t ac_snpB_controls[3] = {0,0,0};
    uint32_t ac_snpB_cases[3] = {0,0,0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; j++) {
            ac_snpA_controls[i] += controls[TO2(i,j)];
            ac_snpA_cases[i] += cases[TO2(i,j)];
            ac_snpB_controls[j] += controls[TO2(i,j)];
            ac_snpB_cases[j] += cases[TO2(i,j)];
        }
    }

    pi_jk[0] = ac_snpA_cases[0]+ac_snpA_controls[0]; //sumsCaseAA1+sumsCtrlAA1;
    pi_jk[2] = ac_snpA_cases[2]+ac_snpA_controls[2]; //sumsCaseaa1+sumsCtrlaa1;
    pi_jk[1] = ac_snpA_cases[1]+ac_snpA_controls[1]; //numInds-pi_jk[0]-pi_jk[2];
    pi_jk[3] = ac_snpB_cases[0]+ac_snpB_controls[0]; //sumsCaseAA2+sumsCtrlAA2;
    pi_jk[5] = ac_snpB_cases[2]+ac_snpB_controls[2]; //sumsCaseaa2+sumsCtrlaa2;
    pi_jk[4] = ac_snpB_cases[1]+ac_snpB_controls[1]; //numInds-pi_jk[3]-pi_jk[5];

    for(int i=0; i<6; i++){
        if(!pi_jk[i]){
            pi_jk[i]++;
        }
    }

    GPUEngine::score_type aux1 = cases[1]+cases[4]+cases[7];
    GPUEngine::score_type aux2 = cases[3]+cases[4]+cases[5];

    uint32_t sumsCaseAA1 = ac_snpA_cases[0];
    uint32_t sumsCaseAA2 = ac_snpB_cases[0];
    uint32_t sumsCaseaa1 = ac_snpA_cases[2];
    uint32_t sumsCaseaa2 = ac_snpB_cases[2];

    uint32_t sumsCtrlAA1 = ac_snpA_controls[0];
    uint32_t sumsCtrlAA2 = ac_snpB_controls[0];
    uint32_t sumsCtrlaa1 = ac_snpA_controls[2];
    uint32_t sumsCtrlaa2 = ac_snpB_controls[2];

    // Firstly the cases
    // p_000 = p_00x*p_0x0*p_x00/...
    p[0] = invCases;
    p[0] *= sumsCaseAA1;
    p[0] *= sumsCaseAA2;
    p[0] *= cases[0]+controls[0];
    p[0] /= pi_jk[0]*pi_jk[3];
    // p_001 = p_00x*p_0x1*p_x01/...
    p[1] = invCases;
    p[1] *= sumsCaseAA1;
    p[1] *= aux1;
    p[1] *= cases[1]+controls[1];
    p[1] /= pi_jk[0]*pi_jk[4];
    // p_002 = p_00x*p_0x2*p_x02/...
    p[2] = invCases;
    p[2] *= sumsCaseAA1;
    p[2] *= sumsCaseaa2;
    p[2] *= cases[2]+controls[2];
    p[2] /= pi_jk[0]*pi_jk[5];
    // p_010 = p_01x*p_0x0*p_x10/...
    p[3] = invCases;
    p[3] *= aux2;
    p[3] *= sumsCaseAA2;
    p[3] *= cases[3]+controls[3];
    p[3] /= pi_jk[1]*pi_jk[3];
    // p_011 = p_01x*p_0x1*p_x11/...
    p[4] = invCases;
    p[4] *= aux2;
    p[4] *= aux1;
    p[4] *= cases[4]+controls[4];
    p[4] /= pi_jk[1]*pi_jk[4];
    // p_012 = p_01x*p_0x2*p_x12/...
    p[5] = invCases;
    p[5] *= aux2;
    p[5] *= sumsCaseaa2;
    p[5] *= cases[5]+controls[5];
    p[5] /= pi_jk[1]*pi_jk[5];
    // p_020 = p_02x*p_0x0*p_x20/...
    p[6] = invCases;
    p[6] *= sumsCaseaa1;
    p[6] *= sumsCaseAA2;
    p[6] *= cases[6]+controls[6];
    p[6] /= pi_jk[2]*pi_jk[3];
    // p_021 = p_02x*p_0x1*p_x21/...
    p[7] = invCases;
    p[7] *= sumsCaseaa1;
    p[7] *= aux1;
    p[7] *= cases[7]+controls[7];
    p[7] /= pi_jk[2]*pi_jk[4];
    // p_022 = p_02x*p_0x2*p_x22/...
    p[8] = invCases;
    p[8] *= sumsCaseaa1;
    p[8] *= sumsCaseaa2;
    p[8] *= cases[8]+controls[8];
    p[8] /= pi_jk[2]*pi_jk[5];

    aux1 = controls[1]+controls[4]+controls[7];
    aux2 = controls[3]+controls[4]+controls[5];

    // Then the controls
    // p_100 = p_10x*p_1x0*p_x00/...
    p[9] = invCtrls;
    p[9] *= sumsCtrlAA1;
    p[9] *= sumsCtrlAA2;
    p[9] *= cases[0]+controls[0];
    p[9] /= pi_jk[0]*pi_jk[3];
    // p_101 = p_10x*p_1x1*p_x01/...
    p[10] = invCtrls;
    p[10] *= sumsCtrlAA1;
    p[10] *= aux1;
    p[10] *= cases[1]+controls[1];
    p[10] /= pi_jk[0]*pi_jk[4];
    // p_102 = p_10x*p_1x2*p_x02/...
    p[11] = invCtrls;
    p[11] *= sumsCtrlAA1;
    p[11] *= sumsCtrlaa2;
    p[11] *= cases[2]+controls[2];
    p[11] /= pi_jk[0]*pi_jk[5];
    // p_110 = p_11x*p_1x0*p_x10/...
    p[12] = invCtrls;
    p[12] *= aux2;
    p[12] *= sumsCtrlAA2;
    p[12] *= cases[3]+controls[3];
    p[12] /= pi_jk[1]*pi_jk[3];
    // p_111 = p_11x*p_1x1*p_x11/...
    p[13] = invCtrls;
    p[13] *= aux2;
    p[13] *= aux1;
    p[13] *= cases[4]+controls[4];
    p[13] /= pi_jk[1]*pi_jk[4];
    // p_112 = p_11x*p_1x2*p_x12/...
    p[14] = invCtrls;
    p[14] *= aux2;
    p[14] *= sumsCtrlaa2;
    p[14] *= cases[5]+controls[5];
    p[14] /= pi_jk[1]*pi_jk[5];
    // p_120 = p_12x*p_1x0*p_x20/...
    p[15] = invCtrls;
    p[15] *= sumsCtrlaa1;
    p[15] *= sumsCtrlAA2;
    p[15] *= cases[6]+controls[6];
    p[15] /= pi_jk[2]*pi_jk[3];
    // p_121 = p_12x*p_1x1*p_x21/...
    p[16] = invCtrls;
    p[16] *= sumsCtrlaa1;
    p[16] *= aux1;
    p[16] *= cases[7]+controls[7];
    p[16] /= pi_jk[2]*pi_jk[4];
    // p_122 = p_12x*p_1x2*p_x22/...
    p[17] = invCtrls;
    p[17] *= sumsCtrlaa1;
    p[17] *= sumsCtrlaa2;
    p[17] *= cases[8]+controls[8];
    p[17] /= pi_jk[2]*pi_jk[5];

    // Calculate 1/n
    GPUEngine::score_type n = 0.0;
    for(int i=0; i<18; i++){
        n+=p[i];
    }

    // Save the division as it is more expensive than the multiplication
    n = 1/n;

    GPUEngine::score_type value = 0.0;
    for(int i=0; i<9; i++){
        if(cases[i]){
            p[i]*=n;
            aux1 = cases[i]/(p[i]*numInds);
            aux1 = log(aux1);
            value +=  cases[i]*aux1;
        }
    }
    for(int i=0; i<9; i++){
        if(controls[i]){
            p[i+9]*=n;
            aux1 = controls[i]/(p[i+9]*numInds);
            aux1 = log(aux1);
            value +=  controls[i]*aux1;
        }
    }

    return value;
}

// log-linear filter for 2way BOOST and 2way log-linear kernel

__device__ __host__ void filterLog(const uint32_t* cases, const uint32_t* controls, GPUEngine::score_type *result_score, GPUEngine::score_type *result_error, int *result_iterations){

	// Start detection of mu
    GPUEngine::score_type mutmp[18];
    GPUEngine::score_type mu0tmp[18];
    for(int i=0; i<18; i++){
        mutmp[i] = 1.0;
        mu0tmp[i] = 0.0;
    }

    GPUEngine::score_type mu_ij[9];
    GPUEngine::score_type mu_ik[6];
    GPUEngine::score_type mu_jk[6];

    uint32_t values[18];
    for(int i = 0; i < 18; i++)
        if(i < 9)
            values[i] = cases[i];
        else
            values[i] = controls[i-9];

    // We need a variation of the order of values
    uint32_t aux;
    aux = values[1];
    values[1] = values[9];
    values[9] = values[13];
    values[13] = values[15];
    values[15] = values[16];
    values[16] = values[8];
    values[8] = values[4];
    values[4] = values[2];
    values[2] = aux;
    aux = values[3];
    values[3] = values[10];
    values[10] = values[5];
    values[5] = values[11];
    values[11] = values[14];
    values[14] = values[7];
    values[7] = values[12];
    values[12] = values[6];
    values[6] = aux;

    //Iterative Proportional Fitting homogeneous model (see section 8.7.2 of Categorical Data Analysis (2nd Ed.))
    GPUEngine::score_type muError = 18.0;
    GPUEngine::score_type aux1;

    int iter = 0;
    for (iter = 0; iter <= maxLoglinIter && (muError > 0.001); iter++) {
        memcpy(mu0tmp, mutmp, 18*sizeof(GPUEngine::score_type)); //mu0tmp = mutmp;

        // mu_ij
        for(int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                mu_ij[3*i+j] = mutmp[6*i+2*j] + mutmp[6*i+2*j+1];
            }
        }

        //mu_ijk = mu_ijk*n_ij/mu_ij
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<2; k++){
                    if(mu_ij[3*i+j] > 0){
                        mutmp[6*i+2*j+k] /= mu_ij[3*i+j];
                        mutmp[6*i+2*j+k] *= (values[6*i+2*j]+values[6*i+2*j+1]);
                    }
                    else{
                        mutmp[6*i+2*j+k] = 0;
                    }
                }
            }
        }
        // mu_ik
        for(int i=0; i<3; i++){
            for(int k=0; k<2; k++){
                mu_ik[2*i+k] = mutmp[6*i+k]+mutmp[6*i+2+k]+mutmp[6*i+4+k];
            }
        }
        //mu_ijk = mu_ijk*n_ik/mu_ik
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<2; k++){
                    if(mu_ik[2*i+k] > 0){
                        mutmp[6*i+2*j+k] /= mu_ik[2*i+k];
                        mutmp[6*i+2*j+k] *= (values[k+6*i]+values[k+6*i+2]+values[k+6*i+4]);
                    }
                    else{
                        mutmp[6*i+2*j+k] = 0;
                    }
                }
            }
        }
        // mu_jk
        for(int j=0; j<3; j++){
            for(int k=0; k<2; k++){
                mu_jk[2*j+k] = mutmp[2*j+k]+mutmp[6+2*j+k]+mutmp[12+2*j+k];
            }
        }
        //mu_ijk = mu_ijk*n_jk/mu_jk
        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){
                for(int k=0; k<2; k++){
                    if(mu_jk[2*j+k] > 0){
                        mutmp[6*i+2*j+k] /= mu_jk[2*j+k];
                        mutmp[6*i+2*j+k] *= (values[k+2*j]+values[k+2*j+6]+values[k+2*j+12]);
                    }
                    else{
                        mutmp[6*i+2*j+k] = 0;
                    }
                }
            }
        }

        //calculate Error
        muError = 0.0;
        for(int i=0; i<18; i++){
            aux1 = mutmp[i]-mu0tmp[i];
            if(aux1 < 0){
                aux1 *= -1;
            }
            muError += aux1;
        }
    }

    // Now mu is calculated so we can compute the value
    GPUEngine::score_type value = 0.0;
    for(int i=0; i<18; i++){
        if(values[i]){
            aux1 = values[i]/mutmp[i];
            aux1 = log(aux1);
            value += values[i]*aux1;
        }
    }

    *result_score = 2*value;
    *result_error = muError;
    *result_iterations = iter;
}


// 2WAY BOOST KERNEL

__device__ __host__ void BOOSTKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 9;

    // Filter
    GPUEngine::score_type valueFilterKSA;
    GPUEngine::score_type valueFilter;
    GPUEngine::score_type valueError;
    int    valueIter;
    GPUEngine::score_type invCases = 1.0 / numCases;
    GPUEngine::score_type invCtrls = 1.0 / numControls;

    // KSA
    valueFilterKSA = filterKSA(cases, controls, numCases, numControls, invCases, invCtrls);

    if(valueFilterKSA < cudaKSAThres) {    // Discarded by the KSA filter

        valueFilter = 0.0;
        valueError = 0.0;
        valueIter = 0;

    } else {

        // log-linear test
        filterLog(cases, controls, &valueFilter, &valueError, &valueIter);

        // unnecessary
//        if(2*valueFilter < _cudaLogThres){    // Discarded by the Log Linear filter
//            valueFilter = 0.0;
//        }
    }

    // write back results
    result.setScore(nextResultIndex++, valueFilter);
    result.setScore(nextResultIndex++, valueError);
    result.setScore(nextResultIndex++, (valueFilter > 0 ? approx_p_4df(valueFilter) : 1.0)); // approximate p-value assuming a chi-squared distribution with 4 df

    if (decomp) {
        result.setScore(nextResultIndex++, (GPUEngine::score_type) valueIter);
        result.setScore(nextResultIndex++, valueFilterKSA);
    }


}

// END OF 2WAY BOOST KERNEL


// 2WAY LOG-LINEAR KERNEL

__device__ __host__ void LogLinearKernel2WayCore(const uint32_t *ctable,
                                        uint32_t numCases __attribute__ ((unused)),
                                        uint32_t numControls __attribute__ ((unused)),
                                        DeviceResult<> result,
                                        int &nextResultIndex,
                                        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 9;

    GPUEngine::score_type valueFilter;
    GPUEngine::score_type valueError;
    int    valueIter;
    filterLog(cases, controls, &valueFilter, &valueError, &valueIter);

    // write back results
    result.setScore(nextResultIndex++, valueFilter);
    result.setScore(nextResultIndex++, valueError);
    result.setScore(nextResultIndex++, approx_p_4df(valueFilter)); // approximate p-value assuming a chi-squared distribution with 4 df

    if (decomp) {
        result.setScore(nextResultIndex++, (GPUEngine::score_type) valueIter);
    }

}

// END OF 2WAY LOG-LINEAR KERNEL


// helper functions for logistic regression kernels

__device__ __host__ void fill_zero(double *array, int cnt) {
    for (int i = 0; i < cnt; i++)
        *array++ = 0.0;
}

__device__ __host__ void fill_zero(float *array, int cnt) {
    for (int i = 0; i < cnt; i++)
        *array++ = 0.0;
}

__device__ __host__ void compute_p_and_v_2way(const GPUEngine::score_type* beta, const uint32_t* nijk, GPUEngine::score_type* pijk, GPUEngine::score_type* vij){
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // multiplicative (PLINK)
            GPUEngine::score_type tmp = beta[0] + i*beta[1] + j*beta[2] + i*j*beta[3]; // \sum_j beta_j
            tmp = 1.0/(1.0 + exp(-tmp)); // p_ij
            pijk[TO2(i,j)] = tmp; // controls
            pijk[TO2(i,j)+9] = tmp - 1.0; // cases
            vij[TO2(i,j)] = tmp*(1.0-tmp)*(nijk[TO2(i,j)]+nijk[TO2(i,j)+9]);
        }
    }
}

__device__ __host__ void compute_p_and_v_3way(const GPUEngine::score_type* beta, const uint32_t* nijkl, GPUEngine::score_type* pijkl, GPUEngine::score_type* vijk){
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                // \sum_j beta_j
                GPUEngine::score_type tmp = beta[0] + i*beta[1] + j*beta[2] + k*beta[3]
                           + i*j*beta[4] + i*k*beta[5] + j*k*beta[6] + i*j*k*beta[7];
                tmp = 1.0/(1.0 + exp(-tmp)); // p_ij
                pijkl[TO3(i,j,k)] = tmp; // controls
                pijkl[TO3(i,j,k)+27] = tmp - 1.0; // cases
                vijk[TO3(i,j,k)] = tmp*(1.0-tmp)*(nijkl[TO3(i,j,k)]+nijkl[TO3(i,j,k)+27]);
            }
        }
    }
}

__device__ __host__ void compute_gradient_2way(const uint32_t* nijk, const GPUEngine::score_type* pijk, GPUEngine::score_type* grad){
    const uint32_t *nctrlptr = nijk;
    const uint32_t *ncaseptr = &nijk[9];
    const GPUEngine::score_type *pctrlptr = pijk;
    const GPUEngine::score_type *pcaseptr = &pijk[9];
    fill_zero(grad, 4);
    for (uint32_t i=0; i<3; i++)
        for (uint32_t j=0; j<3; j++) {
            GPUEngine::score_type tmp = (*pctrlptr) * (*nctrlptr) + (*pcaseptr) * (*ncaseptr);
            grad[0] += tmp;
            grad[1] += i*tmp;
            grad[2] += j*tmp;
            grad[3] += i*j*tmp;
            pctrlptr++;
            pcaseptr++;
            nctrlptr++;
            ncaseptr++;
        }
}

__device__ __host__ void compute_gradient_3way(const uint32_t* nijkl, const GPUEngine::score_type* pijkl, GPUEngine::score_type* grad){
    const uint32_t *nctrlptr = nijkl;
    const uint32_t *ncaseptr = &nijkl[27];
    const GPUEngine::score_type *pctrlptr = pijkl;
    const GPUEngine::score_type *pcaseptr = &pijkl[27];
    fill_zero(grad, 8);
    for (uint32_t i=0; i<3; i++)
        for (uint32_t j=0; j<3; j++)
            for (uint32_t k=0; k<3; k++) {
                GPUEngine::score_type tmp = (*pctrlptr) * (*nctrlptr) + (*pcaseptr) * (*ncaseptr);
                grad[0] += tmp;
                grad[1] += i*tmp;
                grad[2] += j*tmp;
                grad[3] += k*tmp;
                grad[4] += i*j*tmp;
                grad[5] += i*k*tmp;
                grad[6] += j*k*tmp;
                grad[7] += i*j*k*tmp;
                pctrlptr++;
                pcaseptr++;
                nctrlptr++;
                ncaseptr++;
            }
}

__device__ __host__ void compute_hessian_2way(const GPUEngine::score_type* vij, GPUEngine::score_type* hh){
    const GPUEngine::score_type *vptr = vij;
    fill_zero(hh, 16);
    for (uint32_t i=0; i<3; i++)
        for (uint32_t j=0; j<3; j++) {
            // h00
            hh[0] += *vptr;
            // h10
            hh[4] += *vptr * i;
            // h11
            hh[5] += *vptr * (i * i);
            // h20
            hh[8] += *vptr * j;
            // h21
            hh[9] += *vptr * (i * j);
            // h22
            hh[10] += *vptr * (j * j);
            // h30 -> == h21
            //hh[12] += *vptr * (i * j);
            // h31
            hh[13] += *vptr * (i * i * j);
            // h32
            hh[14] += *vptr * (i * j * j);
            // h33
            hh[15] += *vptr * (i * i * j * j);
            vptr++;
        }
    hh[12] = hh[9];
}

__device__ __host__ void compute_hessian_3way(const GPUEngine::score_type* vijk, GPUEngine::score_type* hh){
    const GPUEngine::score_type *vptr = vijk;
    fill_zero(hh, 64);
    for (uint32_t i=0; i<3; i++)
        for (uint32_t j=0; j<3; j++)
            for (uint32_t k=0; k<3; k++) {
                // h00
                hh[0] += *vptr;
                // h10
                hh[8] += *vptr * i;
                // h11
                hh[9] += *vptr * (i * i);
                // h20
                hh[16] += *vptr * j;
                // h21
                hh[17] += *vptr * (i * j);
                // h22
                hh[18] += *vptr * (j * j);
                // h30
                hh[24] += *vptr * k;
                // h31
                hh[25] += *vptr * (i * k);
                // h32
                hh[26] += *vptr * (j * k);
                // h33
                hh[27] += *vptr * (k * k);
                // h40 == h21
                //hh[32] += *vptr * (i * j);
                // h41
                hh[33] += *vptr * (i * i * j);
                // h42
                hh[34] += *vptr * (i * j * j);
                // h43
                hh[35] += *vptr * (i * j * k);
                // h44
                hh[36] += *vptr * (i * i * j * j);
                // h50 == h31
                //hh[40] += *vptr * (i * k);
                // h51
                hh[41] += *vptr * (i * i * k);
                // h52 == h43
                //hh[42] += *vptr * (i * j * k);
                // h53
                hh[43] += *vptr * (i * k * k);
                // h54
                hh[44] += *vptr * (i * i * j * k);
                // h55
                hh[45] += *vptr * (i * i * k * k);
                // h60 == h32
                //hh[48] += *vptr * (j * k);
                // h61 == h52 == h43
                //hh[49] += *vptr * (i * j * k);
                // h62
                hh[50] += *vptr * (j * j * k);
                // h63
                hh[51] += *vptr * (j * k * k);
                // h64
                hh[52] += *vptr * (i * j * j * k);
                // h65
                hh[53] += *vptr * (i * j * k * k);
                // h66
                hh[54] += *vptr * (j * j * k * k);
                // h70 == h61 == h52 == h43
                //hh[56] += *vptr * (i * j * k);
                // h71 == h54
                //hh[57] += *vptr * (i * i * j * k);
                // h72 == h64
                //hh[58] += *vptr * (i * j * j * k);
                // h73 == h65
                //hh[59] += *vptr * (i * j * k * k);
                // h74
                hh[60] += *vptr * (i * i * j * j * k);
                // h75
                hh[61] += *vptr * (i * i * j * k * k);
                // h76
                hh[62] += *vptr * (i * j * j * k * k);
                // h77
                hh[63] += *vptr * (i * i * j * j * k * k);
                vptr++;
            }
    hh[32] = hh[17];
    hh[40] = hh[25];
    hh[42] = hh[35];
    hh[48] = hh[26];
    hh[49] = hh[35];
    hh[56] = hh[35];
    hh[57] = hh[44];
    hh[58] = hh[52];
    hh[59] = hh[53];
}

__device__ __host__ void cholesky_decomposition(const GPUEngine::score_type* aa, GPUEngine::score_type* ll, uint32_t param_cnt) {
  GPUEngine::score_type* ll_ptr;
  GPUEngine::score_type* ll_ptr2;
  GPUEngine::score_type fxx;
  GPUEngine::score_type fyy;
  uint32_t param_cntp1 = param_cnt+1;
  for (uint32_t row_idx = 0; row_idx < param_cnt; row_idx++) {
    fxx = aa[row_idx * param_cntp1]; // diagonal element
    ll_ptr = &(ll[row_idx * param_cnt]);
    for (uint32_t col_idx = 0; col_idx < row_idx; col_idx++) {
      fyy = (*ll_ptr++);
      fxx -= fyy * fyy;
    }
    if (fxx >= 0.0) {
      fyy = sqrt(fxx);
    } else {
      fyy = 1e-6; // aha! so this is the solution for any sqrt of a negative number...
    }
    ll[row_idx * param_cntp1] = fyy; // diagonal element
    fyy = 1.0 / fyy;
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 < param_cnt; row_idx2++) {
      fxx = aa[row_idx2 * param_cnt + row_idx];
      ll_ptr = &(ll[row_idx * param_cnt]);
      ll_ptr2 = &(ll[row_idx2 * param_cnt]);
      for (uint32_t col_idx = 0; col_idx < row_idx; col_idx++) {
        fxx -= (*ll_ptr++) * (*ll_ptr2++);
      }
      ll[row_idx2 * param_cnt + row_idx] = fxx * fyy;
    }
  }
}

__device__ __host__ void cholesky_decomposition_float(const float* aa, float* ll, uint32_t param_cnt) {
  float* ll_ptr;
  float* ll_ptr2;
  float fxx;
  float fyy;
  uint32_t param_cntp1 = param_cnt+1;
  for (uint32_t row_idx = 0; row_idx < param_cnt; row_idx++) {
    fxx = aa[row_idx * param_cntp1]; // diagonal element
    ll_ptr = &(ll[row_idx * param_cnt]);
    for (uint32_t col_idx = 0; col_idx < row_idx; col_idx++) {
      fyy = (*ll_ptr++);
      fxx -= fyy * fyy;
    }
    if (fxx >= 0.0) {
      fyy = sqrtf(fxx);
    } else {
      fyy = 1e-6; // aha! so this is the solution for any sqrt of a negative number...
    }
    ll[row_idx * param_cntp1] = fyy; // diagonal element
    fyy = 1.0 / fyy;
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 < param_cnt; row_idx2++) {
      fxx = aa[row_idx2 * param_cnt + row_idx];
      ll_ptr = &(ll[row_idx * param_cnt]);
      ll_ptr2 = &(ll[row_idx2 * param_cnt]);
      for (uint32_t col_idx = 0; col_idx < row_idx; col_idx++) {
          fxx -= (*ll_ptr++) * (*ll_ptr2++);
      }
      ll[row_idx2 * param_cnt + row_idx] = fxx * fyy;
    }
  }
}

__device__ __host__ void solve_linear_system(const GPUEngine::score_type* ll, const GPUEngine::score_type* yy, GPUEngine::score_type* xx, uint32_t param_cnt) {
  const GPUEngine::score_type* ll_ptr;
  GPUEngine::score_type* xx_ptr;
  uint32_t row_idx;
  uint32_t col_idx;
  GPUEngine::score_type fxx;
  for (row_idx = 0; row_idx < param_cnt; row_idx++) {
    fxx = yy[row_idx];
    ll_ptr = &(ll[row_idx * param_cnt]);
    xx_ptr = xx;
    for (col_idx = 0; col_idx < row_idx; col_idx++) {
      fxx -= (*ll_ptr++) * (*xx_ptr++);
    }
    *xx_ptr = fxx / (*ll_ptr);
  }
  for (col_idx = param_cnt; col_idx;) {
    fxx = xx[--col_idx];
    xx_ptr = &(xx[param_cnt-1]);
    for (row_idx = param_cnt-1; row_idx > col_idx; row_idx--) {
      fxx -= ll[row_idx * param_cnt + col_idx] * (*xx_ptr--);
    }
    *xx_ptr = fxx / ll[row_idx * (param_cnt+1)];
  }
}


__device__ __host__ void solve_linear_system_float(const float* ll, const float* yy, float* xx, uint32_t param_cnt) {
  const float* ll_ptr;
  float* xx_ptr;
  uint32_t row_idx;
  uint32_t col_idx;
  float fxx;
  for (row_idx = 0; row_idx < param_cnt; row_idx++) {
    fxx = yy[row_idx];
    ll_ptr = &(ll[row_idx * param_cnt]);
    xx_ptr = xx;
    for (col_idx = 0; col_idx < row_idx; col_idx++) {
      fxx -= (*ll_ptr++) * (*xx_ptr++);
    }
    *xx_ptr = fxx / (*ll_ptr);
  }
  for (col_idx = param_cnt; col_idx;) {
    fxx = xx[--col_idx];
    xx_ptr = &(xx[param_cnt-1]);
    for (row_idx = param_cnt-1; row_idx > col_idx; row_idx--) {
      fxx -= ll[row_idx * param_cnt + col_idx] * (*xx_ptr--);
    }
    *xx_ptr = fxx / ll[row_idx * (param_cnt+1)];
  }
}

// 2WAY (PLINK) LOGISTIC REGRESSION KERNEL

__device__ __host__ void LogisticRegressionKernel2WayCore(const uint32_t *ctable,
                                        uint32_t numCases __attribute__ ((unused)),
                                        uint32_t numControls __attribute__ ((unused)),
                                        DeviceResult<> result,
                                        int &nextResultIndex,
                                        bool decomp) {

    GPUEngine::score_type beta[4];
    GPUEngine::score_type dbeta[4];
    GPUEngine::score_type pijk[18];
    GPUEngine::score_type vij[9];
    GPUEngine::score_type grad[4];
    GPUEngine::score_type hh[16];
    GPUEngine::score_type ll[16];

    GPUEngine::score_type delta = 0.0;
    GPUEngine::score_type min_delta = 1e9; // according to PLINK...

    int returncode = 0;

    uint32_t iteration = 0;
    fill_zero(beta, 4); // starting point
    fill_zero(hh, 16);
    do {
        iteration++;

        compute_p_and_v_2way(beta, ctable, pijk, vij);
        compute_gradient_2way(ctable, pijk, grad);
        compute_hessian_2way(vij, hh);
        cholesky_decomposition(hh, ll, 4);
        solve_linear_system(ll, grad, dbeta, 4);

        delta = 0.0;
        for (int i = 0; i < 4; i++) {
          delta += fabs(dbeta[i]);
          beta[i] -= dbeta[i];
        }
        if (delta < min_delta) {
          min_delta = delta;
        }
        if (::isnan(delta)) { // test if NaN
            returncode = 1;
            break;
        }
        if (iteration > 4) {
          if ((delta > 20.0) && (delta > 2 * min_delta)) {
              returncode = 2;
              break; // convergence failure
          }
          if ((iteration >= 8) && fabs(1.0 - delta) < 1e-3) {
              returncode = 3;
              break; // convergence failure
          }
        }
    } while(iteration < 16 && delta >= 1e-4); // according to PLINK

    GPUEngine::score_type zsq = beta[3] * beta[3] * ll[15] * ll[15]; // \beta_3^2 / (1/l_33^2)
    GPUEngine::score_type oddr = exp(beta[3]);

    // write back results
    result.setScore(nextResultIndex++, returncode == 0? zsq : 0);
    result.setScore(nextResultIndex++, returncode == 0? oddr : 0);
    result.setScore(nextResultIndex++, returncode == 0? approx_p_1df(zsq) : 1.0); // calculate approximate p-value (assuming chi-square with 1 df)
    result.setScore(nextResultIndex++, returncode == 0? beta[3] : 0);
    result.setScore(nextResultIndex++, returncode == 0? fabs(1/ll[15]) : 0);

    if (decomp) {
//        assert(result.getScoreFields() >= 8);
        result.setScore(nextResultIndex++, (GPUEngine::score_type) iteration);
        result.setScore(nextResultIndex++, delta);
        result.setScore(nextResultIndex++, min_delta);
        result.setScore(nextResultIndex++, (GPUEngine::score_type) returncode);
    }

}

// END OF 2WAY (PLINK) LOGISTIC REGRESSION KERNEL


// 3WAY LOGISTIC REGRESSION KERNEL

__device__ __host__ void LogisticRegressionKernel3WayCore(const uint32_t *ctable,
                                            uint32_t numCases __attribute__ ((unused)),
                                            uint32_t numControls __attribute__ ((unused)),
                                            DeviceResult<> result,
                                            int &nextResultIndex,
                                            bool decomp) {

    GPUEngine::score_type beta[8];
    GPUEngine::score_type dbeta[8];
    GPUEngine::score_type pijkl[54];
    GPUEngine::score_type vijk[27];
    GPUEngine::score_type grad[8];
    GPUEngine::score_type hh[64];
    GPUEngine::score_type ll[64];

    GPUEngine::score_type delta = 0.0;
    GPUEngine::score_type min_delta = 1e9; // according to PLINK...

    int returncode = 0;

    uint32_t iteration = 0;
    fill_zero(beta, 8); // starting point
    do {
        iteration++;

        compute_p_and_v_3way(beta, ctable, pijkl, vijk);
        compute_gradient_3way(ctable, pijkl, grad);
        compute_hessian_3way(vijk, hh);
        cholesky_decomposition(hh, ll, 8);
        solve_linear_system(ll, grad, dbeta, 8);

        delta = 0.0;
        for (int i = 0; i < 8; i++) {
          delta += fabs(dbeta[i]);
          beta[i] -= dbeta[i];
        }
        if (delta < min_delta) {
          min_delta = delta;
        }
        if (::isnan(delta)) { // test if NaN
            returncode = 1;
            break;
        }
        if (iteration > 4) {
          if ((delta > 20.0) && (delta > 2 * min_delta)) {
              returncode = 2;
              break; // convergence failure
          }
          if ((iteration >= 8) && fabs(1.0 - delta) < 1e-3) {
              returncode = 3;
              break; // convergence failure
          }
        }
    } while(iteration < 16 && delta >= 1e-4); // according to PLINK

    GPUEngine::score_type zsq = beta[7] * beta[7] * ll[63] * ll[63]; // \beta_7^2 / (1/l_77^2)
    GPUEngine::score_type oddr = exp(beta[7]);

    // write back results
    result.setScore(nextResultIndex++, returncode == 0? zsq : 0);
    result.setScore(nextResultIndex++, returncode == 0? oddr : 0);
    result.setScore(nextResultIndex++, returncode == 0? approx_p_1df(zsq) : 1.0); // calculate approximate p-value (assuming chi-square with 1 df)
    result.setScore(nextResultIndex++, returncode == 0? beta[7] : 0);
    result.setScore(nextResultIndex++, returncode == 0? fabs(1/ll[63]) : 0);

    if (decomp) {
//        assert(result.getScoreFields() >= 8);
        result.setScore(nextResultIndex++, (GPUEngine::score_type) iteration);
        result.setScore(nextResultIndex++, delta);
        result.setScore(nextResultIndex++, min_delta);
        result.setScore(nextResultIndex++, (GPUEngine::score_type) returncode);
    }
}

// END OF 3WAY LOGISTIC REGRESSION KERNEL


// 2WAY LD KERNEL

__device__ __host__ void LDKernel2WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 9;

    double rsq_high, rsq_low;
    calcLD_2way(cases, controls, numCases, numControls, &rsq_high, &rsq_low); //, id[SNP_A], id[SNP_B]);
    result.setScore(nextResultIndex++, static_cast<GPUEngine::score_type>(rsq_high));
    result.setScore(nextResultIndex++, static_cast<GPUEngine::score_type>(rsq_low));

    if (decomp) { // unused
    }
}

// END OF 2WAY LD KERNEL

// 3WAY LD KERNEL

__device__ __host__ void LDKernel3WayCore(const uint32_t *ctable,
        uint32_t numCases,
        uint32_t numControls,
        DeviceResult<> result,
        int &nextResultIndex,
        bool decomp) {

    const uint32_t *controls = ctable;
    const uint32_t *cases = ctable + 27;

    double ldAB, ldAC, ldBC;
    calcPairwiseLD_3way(cases, controls, numCases, numControls, &ldAB, &ldAC, &ldBC); //, id[SNP_A], id[SNP_B], id[SNP_C]);
    result.setScore(nextResultIndex++, ldAB);
    result.setScore(nextResultIndex++, ldAC);
    result.setScore(nextResultIndex++, ldBC);

    if (decomp) { // unused
    }
}

// END OF 3WAY LD KERNEL

#undef TO2
#undef TO3
#undef TO2A
#undef TO3A

