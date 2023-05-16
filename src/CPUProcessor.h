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

#ifndef CPUPROCESSOR_H_
#define CPUPROCESSOR_H_

#include "ResultView.h"
#include "ResultHandler.h"
#include "Processor.h"

class CPUProcessor : public Processor {
public:
    CPUProcessor() = delete;
    CPUProcessor(const SNPDB &db, const vector<Method> &methods, Progress &progress, ResultView<> &view, ResultHandler<> &rh, bool debug = false);
    ~CPUProcessor();

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
    ResultView<> &view;
    ResultHandler<> &rh;

    size_t resultsize;
    vector<char*> resultbufs;
};

#endif /* CPUPROCESSOR_H_ */
