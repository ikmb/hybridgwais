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

#ifndef FPGAHANDLER_H
#define FPGAHANDLER_H

#include <vector>
#include <list>
#include <atomic>
#include <memory>
#include <chrono>
#include <iostream>
#include <utility>
#include <time.h>
#include <mutex>
#include <condition_variable>
#include <tbb/concurrent_queue.h>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/Buffer.h"
#include "hybridsys/ThreadUtils.h"

#include "SNPDB.h"
#include "FPGAConfigurationGWAIS.h"
#include "Method.h"
#include "ThreadPool.h"
#include "Processor.h"
#include "Progress.h"

using namespace std;

class FPGAHandler : public Processor {
public:

    FPGAHandler(
            hybridsys::Hybridsys &hybridsys,
            chrono::milliseconds timeout,
            hybridsys::BufferFactory<hybridsys::FPGABuffer> &factory,
            const SNPDB &snpdb,
            const vector<Method> &methods,
            Progress &progress,
            const bool &term_request,
            bool debug);

    void createTables();

    void fetchTables(tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &inqueue,
                 tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &outqueue,
                 atomic_bool& termination_request,
                 int threadIndex
                 );

    bool isBlockFinished(unsigned fpgaidx) {
    	bool ret = true;
    	for (size_t br : buffersRemaining[fpgaidx])
    		ret = ret && br == 0;
    	return ret;
    }

    bool isFinished() {
        return finished_flag;
    }

//    // shows the progress of the FPGA buffers between 0 (not started) and 1 (finished)
//    double progress() {
//    	size_t bufrem = 0;
//    	for (size_t br : buffersRemaining)
//    		bufrem += br;
//    	return 1.0 - ((double)bufrem / totalTransferBuffers);
//    }

//    // the fraction of processed buffer count to totalTransferBuffers is the progress
//    size_t getProcessedBuffersCnt() {
//        size_t bufrem = 0;
//        for (size_t br : buffersRemaining)
//            bufrem += br;
//        return totalTransferBuffers - bufrem;
//    }

    // returns true if a timeout has occurred on one channel
    bool timeoutOccurred() {
    	return timeout_flag;
    }

    static constexpr int num_constants = 12;

    static constexpr size_t GPUONLY_IDSIZE_2WAY = 8;
    static constexpr size_t GPUONLY_IDSIZE_3WAY = 12;

    class Block {
    public:
        Block() {}
        Block(const Block &b) = default;
        vector<snprange> sets; // the SNP intervals that are subsequently sent to the FPGA: no overlaps allowed!!!
        size_t init_last = 0; // FPGA stops after initializing and processing this SNP (FPGA local index)
        size_t stream_start = 0; // after initializing the PE chain streaming starts with this SNP (if it was not part of the initialization) (FPGA local index)
        size_t snpcount = 0; // how many SNPs are in the sets?
        vector<size_t> idmap; // mapping of the FPGA local SNPs to the local SNP indices in the SNPDB
        bool swapped = false; // result IDs outside the "triangle" region should be swapped
        size_t excluderange = 0; // for filtering SNP pairs that are in the user defined exclude range during result processing
        vector<size_t> buffercount; // pre-calculated number of result buffers for each PE chain resulting from the specified sets
        vector<size_t> residualresults; // pre-calculated number of results in the last result buffer for each PE chain
        size_t resultcount_filtered_per_buffer = 0; // pre-calculate number of filtered(!) results from this block averaged per buffer (used for the progress)
    };

protected:
    // virtual methods from parent class to process a combination of sets
    void processSetPair(const snprange &setA_, const snprange &setB_, size_t excluderange_) override;
    void processSetTriple(const snprange &setA_, const snprange &setB_, const snprange &setC_, size_t excluderange_) override;
    // virtual methods from parent class need to be implemented, but are not required here
    void processPair(size_t snpA __attribute((unused)), size_t snpB __attribute((unused))) override {}
    void processTriple(size_t snpA __attribute((unused)), size_t snpB __attribute((unused)), size_t snpC __attribute((unused))) override {}


private:
    hybridsys::Hybridsys &hybridsys;
    unsigned num_fpgas;
    FPGAConfigurationGWAIS fpgaconf;
    chrono::milliseconds timeout;
    vector<bool> setready_flags; // for each FPGA
    bool finished_flag = false;
    bool timeout_flag = false;
    unsigned order;
    hybridsys::BufferFactory<hybridsys::FPGABuffer> &bufferFactory;

    // for the current pair of sets:
    vector<vector<size_t>> buffersRemaining; // for each FPGA and for each chain
    vector<shared_ptr<Block>> currentBlocks; // one for each FPGA

    // mutexes and cvs
    vector<vector<mutex*>> mtxs; // one for each FPGA and channel
    vector<vector<condition_variable*>> cvs; // one for each FPGA and channel
    mutex progress_mx;

    size_t tablesPerBuffer;

    vector<shared_ptr<Block>> initDistribution(const snprange &setA_, const snprange &setB_, size_t excluderange_);
    void initializeFPGA(unsigned fpgaindex);
    void sendSNPData(unsigned fpgaindex);

    const bool &term_request;

    bool debug;
};
#endif // FPGAHANDLER_H
