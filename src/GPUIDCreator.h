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

#ifndef GPUIDCREATOR_H_
#define GPUIDCREATOR_H_

#include <tbb/concurrent_queue.h>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/Buffer.h"
#include "hybridsys/ThreadUtils.h"

#include "Processor.h"

class GPUIDCreator : public Processor {
public:
    GPUIDCreator() = delete;
    GPUIDCreator(hybridsys::BufferFactory<hybridsys::FPGABuffer> &factory,
            const SNPDB &db, const vector<Method> &methods, Progress &progress, bool debug = false);
    ~GPUIDCreator() {}

    void createIDs(tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &inqueue __attribute__ ((unused)),
             tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &outqueue,
             atomic_bool& termination_request,
             int threadIndex
             );

    bool isFinished() { return isfinished; }

    static constexpr size_t GPUONLY_IDSIZE_2WAY = 8;
    static constexpr size_t GPUONLY_IDSIZE_3WAY = 12;

protected:
    void processSetPair(const snprange &setA_, const snprange &setB_, size_t excluderange_) override {
    	if (SNPDB::isEmpty(setA_) || SNPDB::isEmpty(setB_))
    		return;
        // call parent function
        processSetPairGeneral(setA_, setB_, excluderange_);
    }
    void processSetTriple(const snprange &setA_, const snprange &setB_, const snprange &setC_, size_t excluderange_) override {
    	if (SNPDB::isEmpty(setA_) || SNPDB::isEmpty(setB_) || SNPDB::isEmpty(setC_))
			return;
        // call parent function
        processSetTripleGeneral(setA_, setB_, setC_, excluderange_);
    }
    void processPair(size_t snpA, size_t snpB) override;
    void processTriple(size_t snpA, size_t snpB, size_t snpC) override;

private:

    void checkBuffer();

    hybridsys::BufferFactory<hybridsys::FPGABuffer> &bufferFactory;
    shared_ptr<hybridsys::FPGABuffer> b; // current buffer
    unsigned *currptr;
    size_t items_in_buffer;
    size_t items_per_buffer;

    tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>>* outqueue_ptr;

    bool isfinished;
};

#endif /* GPUIDCREATOR_H_ */
