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

# include "GPUIDCreator.h"

GPUIDCreator::GPUIDCreator(hybridsys::BufferFactory<hybridsys::FPGABuffer> &factory,
        const SNPDB &db_, const vector<Method> &methods_, Progress &progress_, bool debug_)
    : Processor(db_, methods_, 1, progress_, debug_), // using only 1 thread for ID creation as the process functions are not thread safe! This is definitely a TODO as this is the bottleneck for GPU-only processing!
      bufferFactory(factory)
{
    currptr = NULL;
    items_in_buffer = 0;
    size_t fpgatablesize = order == 2 ? GPUIDCreator::GPUONLY_IDSIZE_2WAY : GPUIDCreator::GPUONLY_IDSIZE_3WAY;
    tables_per_buffer = bufferFactory.getBufferSize() / fpgatablesize;
    outqueue_ptr = NULL;
    isfinished = false;
}

void GPUIDCreator::createIDs(tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &inqueue_ __attribute__ ((unused)),
             tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &outqueue,
             atomic_bool& termination_request,
             int threadIndex __attribute__ ((unused))
             ) {
    ThreadUtils::setThreadName("create IDs");

    items_in_buffer = 0;
    b = bufferFactory.get(); // get a transmission buffer, wait until one is ready
    currptr = (unsigned*) (b->getData());
    outqueue_ptr = &outqueue; // using "address of" at the reference makes a pointer

    // call parent function
    process();

    // send the last buffer
    b->setContentLength(items_in_buffer);
    outqueue_ptr->push(b);

    isfinished = true; // used to stop other processes

    if (termination_request) // should never get here as clean terminations are only required for FPGA runs
        cout << "ID creator terminated." << endl;
    else if (debug)
        cout << "ID creator finished." << endl;

}

void GPUIDCreator::processPair(size_t snpA, size_t snpB) {
    // test snpA and snpB
    checkBuffer();
    *currptr++ = snpA;
    *currptr++ = snpB;
    items_in_buffer++;
}

void GPUIDCreator::processTriple(size_t snpA, size_t snpB, size_t snpC) {
    // test snpA and snpB and snpC
    checkBuffer();
    *currptr++ = snpA;
    *currptr++ = snpB;
    *currptr++ = snpC;
    items_in_buffer++;
}

void GPUIDCreator::checkBuffer() {
    if (items_in_buffer == tables_per_buffer) {
        // buffer is full -> send and get a new one
        b->setContentLength(tables_per_buffer);
        outqueue_ptr->push(b); // send
        if (debug) {
            cout << "Sent full buffer. " << tables_per_buffer << endl;
        }
        b = bufferFactory.get(); // wait for a new one
        items_in_buffer = 0;
        currptr = (unsigned*) (b->getData());
    }
}


// TODO Don't forget that due to an exclude range the number of remaining buffers may shrink!!
