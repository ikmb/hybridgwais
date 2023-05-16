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

#include <omp.h>

#include "CPUProcessor.h"
#include "ResultView.h"
#include "GPUKernels.h"


CPUProcessor::CPUProcessor(const SNPDB &db_, const vector<Method> &methods_, Progress& progress_, ResultView<> &view_, ResultHandler<> &rh_, bool debug_)
    : Processor(db_, methods_, omp_get_max_threads(), progress_, debug_), view(view_), rh(rh_) { // TODO set number of threads to user value?
    resultsize = view_.getResultSize();
    resultbufs.reserve(omp_get_max_threads());
    for (int i = 0; i < omp_get_max_threads(); i++)
        resultbufs.push_back(new char[resultsize]);
}

CPUProcessor::~CPUProcessor() {
    for (const auto &r : resultbufs)
        if (r) delete[] r;
}

void CPUProcessor::processPair(size_t snpA, size_t snpB) {
//    // DEBUG
//    cerr << snpA+1 << "\t" << snpB+1 << endl;

    int tn = omp_get_thread_num();

    // test SNP a and SNP b
    ctable2way t;
    uint32_t numcases, numctrls;
	generateContingencyTable2Way(t.data(), &numcases, &numctrls,
			db.data() + snpA * db.getSNPSize(),
			db.data() + snpB * db.getSNPSize(),
			db.getSNPSize(), db.getCaseCountPadded(), snpA, snpB);
//        // DEBUG
//        printf("\nctable: snpsize: %lu casecnt: %lu dbA: %016lx dbB: %016lx\n", db.getSNPSize(), db.getCaseCountPadded(), (uint64_t) (db.data() + snpA * db.getSNPSize()), (uint64_t) (db.data() + snpB * db.getSNPSize()));
//        for (int i = 0; i < 18; i++) {
//            printf("%d\t", t.data()[i]);
//            if (i%3 == 2)
//                printf("\n");
//            if (i==8)
//                printf("\n");
//        }
//        printf("cases: %d  controls: %d\n", numcases, numctrls);
//        printf("ID A: %lu ID B: %lu\n\n", snpA, snpB);
//        // __DEBUG

    char* resultbuf = resultbufs[tn];
    ResultView<> v(view); // local copy
    v.setBuffer(resultbuf, resultsize);

    setID2Way(t.data(), v.getDeviceResult(0));

    int nextResultIndex = 0;
    for (const auto &method : methods) {
        switch (method.getType()) {
            case Method::MI:
                MutualInformationKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            case Method::IG:
                InformationGainKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            case Method::BOOST:
                BOOSTKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            case Method::LOGLIN:
                LogLinearKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            case Method::LOGREG:
				LogisticRegressionKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            case Method::LD:
                LDKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
                break;
            default: // 3way or unknown kernel (should never get here)
                throw runtime_error("Undefined kernel!");
        }
    }

    rh.processResultView(v, tn);

}

void CPUProcessor::processTriple(size_t snpA, size_t snpB, size_t snpC) {
//    // DEBUG
//    cerr << snpA+1 << "\t" << snpB+1 << "\t" << snpC+1 << endl;

    int tn = omp_get_thread_num();

    // test snpA and snpB and snpC
    ctable3way t;
    uint32_t numcases, numctrls;
    generateContingencyTable3Way(t.data(), &numcases, &numctrls,
            db.data() + snpA * db.getSNPSize(),
            db.data() + snpB * db.getSNPSize(),
            db.data() + snpC * db.getSNPSize(),
            db.getSNPSize(), db.getCaseCountPadded(), snpA, snpB, snpC);

    char* resultbuf = resultbufs[tn];
    ResultView<> v(view); // local copy
    v.setBuffer(resultbuf, resultsize);

    setID3Way(t.data(), v.getDeviceResult(0));

    int nextResultIndex = 0;
    for (const auto &method : methods) {
        switch (method.getType()) {
        case Method::MI3:
            MutualInformationKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
            break;
        case Method::IG3:
            InformationGainKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
            break;
        case Method::LOGREG3:
            LogisticRegressionKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
            break;
        case Method::LD3:
            LDKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), nextResultIndex, method.isDetailedComputation());
            break;
        default: // 2way or unknown kernel (should never get here)
            throw runtime_error("Undefined kernel!");
        }
    }

    rh.processResultView(v, tn);

}
