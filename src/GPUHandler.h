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

#ifndef GPUHANDLER_H
#define GPUHANDLER_H

#include <atomic>
#include <array>
#include <cinttypes>
#include <string>
#include <vector>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/Buffer.h"
#include "hybridsys/ThreadUtils.h"

#include "FPGAConfigurationGWAIS.h"
#include "GPUEngine.h"
#include "ThreadPool.h"


using namespace std;

class GPUEngine;
class SNPDB;

class GPUHandler
{
public:

    using score_type = GPUEngine::score_type;
    using id_type = uint32_t;

    GPUHandler(hybridsys::Hybridsys &hybridsys_, const SNPDB &db, hybridsys::BufferFactory<hybridsys::CUDABuffer> &factory, size_t tableSize, size_t tableBufferSize, const vector<Method> &methods_, ResultView<> view, bool useFPGA_, bool debug_)
        : hybridsys(hybridsys_),
          methods(methods_),
          snpdb(db),
          engines(),
          bufferFactory(factory),
          useFPGA(useFPGA_),
		  debug(debug_)
    {
        int num_gpus = hybridsys.getGPUs().size();
        for(int i = 0; i < num_gpus; i++)
            engines.emplace_back(snpdb, methods, hybridsys.getGPU(i).getIndex(), tableSize, tableBufferSize, view, useFPGA, debug);
    }

    void distributeSNPData() {
        for (auto & e : engines)
            e.uploadSNPData();
    }

    template<typename Tin, typename Tout>
    void process(typename ThreadPool<Tin, Tout>::inqueue_type &inqueue,
                 typename ThreadPool<Tin, Tout>::outqueue_type &outqueue,
                 atomic_bool& termination_request,
                 int threadIndex
                 ) {

        ThreadUtils::setThreadName("GPU:" + to_string(engines[threadIndex].getCUDAIndex()));

        engines[threadIndex].initialize();

        while(!termination_request) {
            Tin b;
            try {
                inqueue.pop(b);
            } catch (tbb::user_abort &e) { // is thrown whenever the process is cancelled or regularly finished
                break;
            }

            auto target = bufferFactory.get();

            const unsigned long tablesExpected = b->getContentLength();

            target->setContentLength(tablesExpected);

            // copy metadata (for ID mapping if necessary) from source buffer
            target->setMetadata(b->getMetadata());

            engines[threadIndex].runKernel(b->getData(), target->getData(), tablesExpected);

            outqueue.push(target);
        }

        engines[threadIndex].free();
    }

    void setDumpStream(ostream &dump) {
        for (auto &e : engines)
            e.setDumpStream(dump);
    }

private:
    hybridsys::Hybridsys &hybridsys;
    const vector<Method> &methods;
    const SNPDB &snpdb;
    vector<GPUEngine> engines;
    hybridsys::BufferFactory<hybridsys::CUDABuffer> &bufferFactory;
    bool useFPGA;
    bool debug;
};



#endif // GPUHANDLER_H
