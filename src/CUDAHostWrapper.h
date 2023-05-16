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

#ifndef CUDAHOSTWRAPPER_H
#define CUDAHOSTWRAPPER_H

#ifndef USE_CUDA_GPU

#define __constant__
#define __host__
#define __device__
#define __global__

#define checkCUDAError(x)
#include <cmath>
#include <cstdint>
#include <cstring>

struct dummy_thread_coords {
    uint32_t x;
    uint32_t y;
};

static struct dummy_thread_coords threadIdx {0, 0};
static struct dummy_thread_coords blockDim {0, 0};
static struct dummy_thread_coords blockIdx {0, 0};

#endif

#endif // CUDAHOSTWRAPPER_H
