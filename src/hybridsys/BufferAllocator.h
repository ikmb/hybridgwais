/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
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

#ifndef CUDAALLOCATOR_H
#define CUDAALLOCATOR_H

#include <cassert>
#include <memory>

#ifdef USE_CUDA_GPU
#include <cuda_runtime.h>
#endif

#include <iostream>
#include <cstdio>

extern "C" {
#include <unistd.h>
}

namespace hybridsys {

using namespace std;

template<typename T>
struct CUDAAllocator {
    using value_type = T;

    // default constructors
    CUDAAllocator() noexcept {}
    template <class U> CUDAAllocator(const CUDAAllocator<U> &) {}

   /* static */ T* allocate(size_t n) {
       T *ptr;

       // check overflow of size_t
       size_t total_size;
       total_size = n * sizeof(T);

#ifdef USE_CUDA_GPU
       cudaError_t status = cudaHostAlloc(&ptr, total_size, 0);
       if(status != cudaSuccess)
           throw bad_alloc();
#else
       ptr = (char*) malloc(total_size);
       if(ptr == NULL)
           throw bad_alloc();
#endif

       return ptr;
    }

    static void deallocate(T* ptr, size_t n __attribute__((unused))) {
#ifdef USE_CUDA_GPU
        cudaFreeHost(ptr);
#else
        free(ptr);
#endif
    }
};

template<typename T>
struct PageAlignedAllocator {
    using value_type = T;

    PageAlignedAllocator() noexcept {}
    template<class U> PageAlignedAllocator(const PageAlignedAllocator<U> &) {}

    /* static */ T* allocate(size_t n) {
        T* ptr;

        size_t total_size;
        total_size = n * sizeof(T);

        size_t page_size = sysconf(_SC_PAGESIZE);
        ptr = reinterpret_cast<T*>(aligned_alloc(page_size, total_size));
        if(!ptr) {
            throw bad_alloc();
        };

#ifdef LOCK_FPGA_BUFFERS
        if(mlock(ptr, total_size) == -1) {
            throw bad_alloc();
        }
#endif
        return ptr;

    }

    static void deallocate(T* ptr, size_t n __attribute__((unused))) {
#ifdef LOCK_FPGA_BUFFERS
        munlock(ptr, n);
#endif
        free(ptr);
    }

};

} // namespace

#endif // CUDAALLOCATOR_H
