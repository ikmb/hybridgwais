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

#ifdef USE_CUDA_GPU

#include <sstream>
#include <ostream>
#include <algorithm>

#include "Args.h"

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
extern "C" {
#include <cuda_runtime.h>
}

#include "Method.h"
#include "ResultView.h"
#include "SNPDB.h"
#include "GPUHandler.h"
#include "GPUEngine.h"
#include "GPUKernels.h"

using namespace std;

static cudaDeviceProp gpuProps;

GPUEngine::GPUEngine(const SNPDB &db, const vector<Method> &methods_, int index_, size_t tablesize, size_t tablebuffersize, ResultView<> view, bool useFPGA_, bool debug_)
    : snpdb(db), methods(methods_), order(methods_[0].getOrder()), index(index_), tableSize(tablesize), tableBufferSize(tablebuffersize), resultView(view), useFPGA(useFPGA_), debug(debug_)
{
    const unsigned long max_tables_per_buffer = tableBufferSize / tableSize;
    resultBufferSize = max_tables_per_buffer * view.getResultSize();

    dbOffset = 0;

    devTables = NULL;
    devResults = NULL;
}

void GPUEngine::initialize() {
    cudaSetDevice(index);

    if(gpuProps.multiProcessorCount == 0) {
        checkCUDAError(cudaGetDeviceProperties(&gpuProps, 0))
    }

    if (debug) {
        cout << "GPU Probs: " << index << endl;
        cout << "MPU count: " << gpuProps.multiProcessorCount << endl;
        cout << "Max threads per block: " << gpuProps.maxThreadsPerBlock << endl;
        cout << "Max threads per MPU: " << gpuProps.maxThreadsPerMultiProcessor << endl;
        cout << "Shared memory per block: " << gpuProps.sharedMemPerBlock << endl;
        cout << "Shared memory per MPU: " << gpuProps.sharedMemPerMultiprocessor << endl;
        cout << "Total global memory: " << gpuProps.totalGlobalMem << endl;
        cout << "Mem bus width: " << gpuProps.memoryBusWidth << endl;
    }

    size_t allocsize = tableBufferSize;
    if (!useFPGA) { // only allocate space for the database and copy to GPU if we want to skip the FPGA part
        dbOffset = tableBufferSize;
        size_t mbw = gpuProps.memoryBusWidth/8; // in bytes
        if (dbOffset % mbw) // round up to multiple of mbw
            dbOffset += mbw - (dbOffset % mbw);
        size_t dbbufsize = snpdb.getBufferSize();
        if (dbbufsize % mbw) // round up to multiple of mbw
            dbbufsize += mbw - (dbbufsize % mbw);
        allocsize = dbOffset + dbbufsize; // the DB will be located behind the table memory
    }
    checkCUDAError(cudaMalloc(&devTables, allocsize))
    checkCUDAError(cudaMalloc(&devResults, resultBufferSize))
    // we do not need more memory, so reserve half of the total mem as heap
    checkCUDAError(cudaDeviceSetLimit(cudaLimitMallocHeapSize, gpuProps.totalGlobalMem/2))

    if (!useFPGA)
        uploadSNPData();

    if (debug) {
        int idx = -1;
        cudaGetDevice(&idx);
        cout << "GPU Allocations: idx: " << index << " (" << idx << ")" << endl;
        cout << "tables: @" << hex << (unsigned long)devTables << ", " << dec << tableBufferSize << " bytes" << endl;
        cout << "results: @" << hex << (unsigned long)devResults << ", " << dec << resultBufferSize << " bytes" << endl;
        if (!useFPGA) {
            cout << "DB: @" << hex << (unsigned long)(devTables+dbOffset) << ", " << dec << snpdb.getBufferSize() << " bytes" << endl;
        }
    }

    unsigned long numcases = snpdb.getCaseCountPadded();
    unsigned long numsamples = snpdb.getSampleCountPadded();
    unsigned long snpsize = snpdb.getSNPSize();

    copyConstantsToDevice(methods, numcases, numsamples, snpsize, tableSize, dbOffset, !useFPGA);

    resultView.setBuffer(devResults, resultBufferSize);
}

void GPUEngine::free() {
    cudaFree(devTables);
    cudaFree(devResults);
}

void GPUEngine::uploadSNPData() {
    // SNP DB
    char *devDB = devTables + dbOffset;
    int idx = -1;
    cudaGetDevice(&idx);
    if (debug) {
        cout << "Upload SNP data.\nIDX: " << index << "(" << idx << ")" << endl;
        cout << "DB: @" << hex << (unsigned long)devDB << dec << endl;
        cout << "SNP size: " << snpdb.getSNPSize() << endl;
    }
    checkCUDAError(cudaMemcpy(devDB, snpdb.data(), snpdb.getBufferSize(), cudaMemcpyHostToDevice))

}

void GPUEngine::runKernel(char *source, char *results, size_t tablesExpected) {

    int blockSize = gpuProps.maxThreadsPerBlock/4; // empirically ;-)

    size_t numBlocks = tablesExpected / blockSize;
    if(tablesExpected % blockSize)
        numBlocks++;

//    if (debug) {
//        cout << "Block size: " << blockSize << endl;
//        cout << "Tables exp: " << tablesExpected << endl;
//        cout << "Num blocks: " << numBlocks << endl;
////    size_t heapsize;
////    cudaDeviceGetLimit(&heapsize, cudaLimitMallocHeapSize);
////    cout << "Heap size: " << heapsize << "\nbs*nb*2*snpsize: " << (blockSize * numBlocks * 2 * snpdb.getSNPSize()) << "\nte*2*snpsize: " << (tablesExpected * 2 * snpdb.getSNPSize()) << endl;
//    }

    dim3 grid(numBlocks, 1);
    dim3 blocks(blockSize);

    checkCUDAError(cudaMemcpy(devTables, source, tablesExpected * tableSize, cudaMemcpyHostToDevice))

    cudaDeviceSynchronize();

    int detailedCompFlags = 0;
    for(auto it = methods.crbegin(); it != methods.crend(); it++) { // reverse!
        detailedCompFlags <<= 1;
        if (it->isDetailedComputation())
            detailedCompFlags |= 0x1;
    }

    if (order == 2) {
        Kernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                     tablesExpected,
                                                     resultView,
                                                     detailedCompFlags);
    } else { // order == 3
        Kernel3Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                 tablesExpected,
                                                 resultView,
                                                 detailedCompFlags);
    }

    cudaDeviceSynchronize();
    checkCUDAError(cudaGetLastError())
    cudaDeviceSynchronize();

    checkCUDAError(cudaMemcpy(results, devResults, resultView.getResultCount() * resultView.getResultSize(), cudaMemcpyDeviceToHost))
    cudaDeviceSynchronize();
}

#endif /* USE_CUDA_GPU */
