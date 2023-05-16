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

#ifndef PROCESSOR_H_
#define PROCESSOR_H_

#include <time.h>

#include "SNPDB.h"
#include "Method.h"
#include "Progress.h"

typedef array<uint32_t, 20> ctable2way; // 18 + ID
typedef array<uint32_t, 57> ctable3way; // 54 + ID

class Processor {
public:
    Processor() = delete;
    virtual ~Processor();
    // processes the user sets by calling processSetPair(Triple) accordingly
    void process();

protected:
    Processor(const SNPDB &db, const vector<Method> &methods, int numthreads, Progress &progress, bool debug = false);

    // virtual methods called by process()
    virtual void processSetPair(const snprange &setA_, const snprange &setB_, size_t excluderange) = 0;
    virtual void processSetTriple(const snprange &setA_, const snprange &setB_, const snprange &setC_, size_t excluderange) = 0;

    // these methods can be used to implement processSetPair/Triple() above -> they call processPair/Triple() for each pair/triple:
    void processSetPairGeneral(const snprange &setA_, const snprange &setB_, size_t excluderange, bool swapped = false);
    void processSetTripleGeneral(const snprange &setA_, const snprange &setB_, const snprange &setC_, size_t excluderange);
    // assumes setA to start before setB starts, and setA to end before setB ends
    // swapped: true if setA and setB were swapped to reach the precondition, this implies processPair() will be called with swapped IDs, but only outside an overlapping area.
    void processSetPairSimple(const snprange &setA_, const snprange &setB_, size_t excluderange, bool swapped);
    // parallelization over A
    void processSetPairSimpleParA(const snprange &setA_, const snprange &setB_, size_t excluderange, bool swapped);
    // parallelization over B
    void processSetPairSimpleParB(const snprange &setA_, const snprange &setB_, size_t excluderange, bool swapped);


    // called by processSetPair(Triple)General()
    virtual void processPair(size_t snpA, size_t snpB) = 0;
    virtual void processTriple(size_t snpA, size_t snpB, size_t snpC) = 0;

    const SNPDB &db;
    const vector<Method> &methods;
    unsigned order;

    int numthreads;

    Progress &progress;
    size_t numtests_processed = 0;

    bool debug;
};

#endif /* PROCESSOR_H_ */
