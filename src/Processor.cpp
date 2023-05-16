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
#include <atomic>
#include <assert.h>
#include <mutex>

#include "Progress.h"
#include "StatusFile.h"

#include "Processor.h"

Processor::Processor(const SNPDB &db_, const vector<Method> &methods_, int numthreads_, Progress &progress_, bool debug_)
    : db(db_), methods(methods_), order(methods_[0].getOrder()), numthreads(numthreads_), progress(progress_), debug(debug_) {}

Processor::~Processor(){}

void Processor::process() {

    // compute the number of all tests
    size_t numtests = 0;

    // The SNPDB has already prepared the required intervals for efficient processing.
    // We just need to calculate the required number of tests respective to these intervals before
    // processing the regions itself.

    const SNPDB::RegionInfo &r = db.getRegionInfo();

	// for the different scenarios we precalculate the required combinations of intervals here
	// note that empty intervals still have consistent borders!
	snprange l01 = make_pair(r.l0.first, r.l1.second);
	snprange l02 = make_pair(r.l0.first, r.l2.second);
	snprange l03 = make_pair(r.l0.first, r.l3.second);
	snprange l04 = make_pair(r.l0.first, r.l4.second);
	snprange l0b1 = make_pair(r.l0b.first, r.l1.second);
//	snprange l0b2 = make_pair(r.l0b.first, r.l1.second);
//	snprange l0b3 = make_pair(r.l0b.first, r.l1.second);
//	snprange l0b4 = make_pair(r.l0b.first, r.l1.second);
//	snprange l12 = make_pair(r.l1.first, r.l2.second);
	snprange l13 = make_pair(r.l1.first, r.l3.second);
	snprange l14 = make_pair(r.l1.first, r.l4.second);
	snprange l23 = make_pair(r.l2.first, r.l3.second);
	snprange l24 = make_pair(r.l2.first, r.l4.second);
	snprange l34 = make_pair(r.l3.first, r.l4.second);

    if (order == 2) {

		// four different cases:
		if (r.l0isA) {
			if (r.l4isA) { // case a
				numtests += SNPDB::getPairsCnt(l01, l13);
				numtests += SNPDB::getPairsCnt(l34, l13);
				numtests -= SNPDB::getPairsCnt(r.l1, r.l3); // double calculation if processed naively
			} else { // case b
				numtests += SNPDB::getPairsCnt(l01, l14);
				numtests += SNPDB::getPairsCnt(r.l3, l24);
			}
		} else {
			if (r.l4isA) { // case c
				numtests += SNPDB::getPairsCnt(r.l1, l02);
				numtests += SNPDB::getPairsCnt(l34, l03);
			} else { // case d
				numtests += SNPDB::getPairsCnt(r.l1, l04);
				numtests += SNPDB::getPairsCnt(r.l3, l04);
				numtests -= SNPDB::getPairsCnt(r.l1, r.l3); // double calculation if processed naively
			}
		}

    } else { // order == 3
    	// we handle sets A and B as for 2nd order tests and simply apply the third set C
    	// but we need to pay attention to the special intervals l0b and l4b!
    	// -> This implies some double calculations but we don't care for now... TODO

    	if (r.l0isA) {
			if (r.l4isA) { // case a
//				numtests += SNPDB::getTriplesCnt(l01, l13, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l0, l13, r.lc);
				numtests += SNPDB::getTriplesCnt(l0b1, l13, r.lc);
				numtests += SNPDB::getTriplesCnt(l34, l23, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l4b, l23, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l4, r.l1, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l4b, r.l1, r.lc);
			} else { // case b
//				numtests += SNPDB::getTriplesCnt(l01, l14, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l0, l14, r.lc);
				numtests += SNPDB::getTriplesCnt(l0b1, l14, r.lc);
				numtests += SNPDB::getTriplesCnt(r.l3, l24, r.lc);
				// no need to consider l4b as l4 is B and hence l4b is empty
			}
		} else {
			if (r.l4isA) { // case c
				numtests += SNPDB::getTriplesCnt(r.l1, l02, r.lc); // valid because l0 is B, so l0b is empty
				numtests += SNPDB::getTriplesCnt(l34, l03, r.lc); // valid because l0 is B, so l0b is empty
				numtests += SNPDB::getTriplesCnt(r.l4b, l03, r.lc); // valid because l0 is B, so l0b is empty
			} else { // case d
				numtests += SNPDB::getTriplesCnt(r.l1, l04, r.lc); // valid because l0 is B, so l0b is empty
				numtests += SNPDB::getTriplesCnt(r.l3, l24, r.lc); // valid because l0 is B, so l0b is empty
				numtests += SNPDB::getTriplesCnt(r.l3, r.l0, r.lc);
				// no need to consider l4b as l4 is B and hence l4b is empty
			}
		}
    }

    if (numtests == 0) {
        StatusFile::addError("Test space is empty!");
        exit(EXIT_FAILURE);
    }

    if (debug)
    	cout << " Num tests:\t" << numtests << endl;

    progress.start(numtests);
    numtests_processed = 0;

    // compute all interactions for all combinations of specified sets

	if (order == 2) { // process pairwise

		if (r.l0isA) {
			if (r.l4isA) { // case a
				processSetPair(l01, l13, r.excluderange);
				processSetPair(l34, l23, r.excluderange);
				processSetPair(r.l4, r.l1, r.excluderange);
			} else { // case b
				processSetPair(l01, l14, r.excluderange);
				processSetPair(r.l3, l24, r.excluderange);
			}
		} else {
			if (r.l4isA) { // case c
				processSetPair(r.l1, l02, r.excluderange);
				processSetPair(l34, l03, r.excluderange);
			} else { // case d
				processSetPair(r.l1, l04, r.excluderange);
				processSetPair(r.l3, l24, r.excluderange);
				processSetPair(r.l3, r.l0, r.excluderange);
			}
		}

	} else { // 3rd order
		// we handle sets A and B as for 2nd order tests and simply apply the third set C
		// but we need to pay attention to the special intervals l0b and l4b!
		// -> This implies some double calculations but we don't care for now... TODO

		if (r.l0isA) {
			if (r.l4isA) { // case a
//				processSetTriple(l01, l13, r.lc, r.excluderange);
				processSetTriple(r.l0, l13, r.lc, r.excluderange);
				processSetTriple(l0b1, l13, r.lc, r.excluderange);
				processSetTriple(l34, l23, r.lc, r.excluderange);
				processSetTriple(r.l4b, l23, r.lc, r.excluderange);
				processSetTriple(r.l4, r.l1, r.lc, r.excluderange);
				processSetTriple(r.l4b, r.l1, r.lc, r.excluderange);
			} else { // case b
//				processSetTriple(l01, l14, r.lc, r.excluderange);
				processSetTriple(r.l0, l14, r.lc, r.excluderange);
				processSetTriple(l0b1, l14, r.lc, r.excluderange);
				processSetTriple(r.l3, l24, r.lc, r.excluderange);
				// no need to consider l4b as l4 is B and hence l4b is empty
			}
		} else {
			if (r.l4isA) { // case c
				processSetTriple(r.l1, l02, r.lc, r.excluderange); // valid because l0 is B, so l0b is empty
				processSetTriple(l34, l03, r.lc, r.excluderange); // valid because l0 is B, so l0b is empty
				processSetTriple(r.l4b, l03, r.lc, r.excluderange); // valid because l0 is B, so l0b is empty
			} else { // case d
				processSetTriple(r.l1, l04, r.lc, r.excluderange); // valid because l0 is B, so l0b is empty
				processSetTriple(r.l3, l24, r.lc, r.excluderange); // valid because l0 is B, so l0b is empty
				processSetTriple(r.l3, r.l0, r.lc, r.excluderange);
				// no need to consider l4b as l4 is B and hence l4b is empty
			}
		}

	}

//    // DEBUG
//    cout << " Counted tests:\t" << numtests_processed << endl;
}


void Processor::processSetPairSimple(const snprange &setA, const snprange &setB, size_t excluderange, bool swapped) {

	// this is the precondition!
	assert(setA.first <= setB.first);
	assert(setA.second <= setB.second);

	// dependent on which set is larger, we could choose the parallelization over setA or setB
	if (setA.second - setA.first > setB.second - setB.first)
		processSetPairSimpleParA(setA, setB, excluderange, swapped);
	else
		processSetPairSimpleParB(setA, setB, excluderange, swapped);
}

void Processor::processSetPairSimpleParA(const snprange &setA, const snprange &setB, size_t excluderange, bool swapped) {

	mutex m;

	omp_set_num_threads(numthreads);
#ifndef DEBUG_SINGLETHREAD
#pragma omp parallel for schedule(dynamic)
#endif
	for (size_t a = setA.first; a < setA.second; a++) {
		string chrA;
		size_t bpposA = 0;
		if (excluderange > 0) { // if an exclude range is defined
			chrA = db.getSNPInfo(a).chromosome;
			bpposA = db.getSNPInfo(a).pos_bp;
		}
        size_t bstart = max(a+1, setB.first);
		for (size_t b = bstart; b < setB.second; b++) {
			string chrB;
			size_t bpposB = 0;
			if (excluderange > 0) { // if an exclude range is defined
				chrB = db.getSNPInfo(b).chromosome;
				bpposB = db.getSNPInfo(b).pos_bp;
			}
			// apply excluderange
			if (excluderange == 0 || chrA.compare(chrB) || bpposB - bpposA >= excluderange) {
				if (swapped && (a < setB.first || b >= setA.second)) // only swap, if necessary, outside the overlap (this is only for a convenient representation of the results)
					processPair(b, a);
				else
					processPair(a, b);
//				// DEBUG
//				m.lock();
//				numtests_processed++;
//				m.unlock();
//				// __DEBUG
			}
		}

		// progress (pairs excluded by an excluderange are still counted - which is necessary for the progress indicator!)
		m.lock();
		numtests_processed += setB.second - bstart;
		progress.updateProgress(numtests_processed);
		m.unlock();
	}
}

void Processor::processSetPairSimpleParB(const snprange &setA, const snprange &setB, size_t excluderange, bool swapped) {

//	// DEBUG
//	mutex m;

	for (size_t a = setA.first; a < setA.second; a++) {
		string chrA;
		size_t bpposA = 0;
		if (excluderange > 0) { // if an exclude range is defined
			chrA = db.getSNPInfo(a).chromosome;
			bpposA = db.getSNPInfo(a).pos_bp;
		}
        size_t bstart = max(a+1, setB.first);
		omp_set_num_threads(numthreads);
#ifndef DEBUG_SINGLETHREAD
#pragma omp parallel for schedule(dynamic)
#endif
		for (size_t b = bstart; b < setB.second; b++) {
			string chrB;
			size_t bpposB = 0;
			if (excluderange > 0) { // if an exclude range is defined
				chrB = db.getSNPInfo(b).chromosome;
				bpposB = db.getSNPInfo(b).pos_bp;
			}
			// apply excluderange
			if (excluderange == 0 || chrA.compare(chrB) || bpposB - bpposA >= excluderange) {
				if (swapped && (a < setB.first || b >= setA.second)) // only swap, if necessary, outside the overlap
					processPair(b, a);
				else
					processPair(a, b);
//				// DEBUG
//				m.lock();
//				numtests_processed++;
//				m.unlock();
//				// __DEBUG
			}
		}

		// progress (pairs excluded by an excluderange are still counted - which is necessary for the progress indicator!)

		numtests_processed += setB.second - bstart;
		progress.updateProgress(numtests_processed);
	}
}

void Processor::processSetPairGeneral(const snprange &setA, const snprange &setB, size_t excluderange, bool swapped) {

	// divide up the sets if necessary to get the most efficient processing scheme by using processSetPairSimple()
	if (setA.first <= setB.first) { // first precondition already ok
		if (setA.second <= setB.second) { // second precondition ok
			// nothing to divide or swap
			processSetPairSimple(setA, setB, excluderange, swapped);
		} else { // second precondition not met
			// cut of the front of A and process separately the front of A against B and then B against the remainder of A
			snprange setA1 = make_pair(setA.first, setB.first);
			snprange setA2 = make_pair(setB.first, setA.second);
			if (!SNPDB::isEmpty(setA1))
				processSetPairSimple(setA1, setB, excluderange, swapped);
			if (!SNPDB::isEmpty(setA2))
				processSetPairSimple(setB, setA2, excluderange,	!swapped);
		}
	} else { // setA.first > setB.first: see if we reach the second precondition by simply swapping the sets
		// we could simply achieve this by calling this function recursively with swapped sets
		processSetPairGeneral(setB, setA, excluderange, true);
	}
}

void Processor::processSetTripleGeneral(const snprange &setA, const snprange &setB, const snprange &setC, size_t excluderange) {
//    // DEBUG
//    unordered_set<size_t> alltriples;
//    size_t allcounter = 0;
//    for (size_t a = setA_.first; a < setA_.second; a++) {
//        for (size_t b = setB_.first; b < setB_.second; b++) {
//            for (size_t c = setC_.first; c < setC_.second; c++) {
//                if (a != b && a != c && b != c && alltriples.emplace(xyz(a, b, c)).second)
//                    allcounter++;
//            }
//        }
//    }
//
//    unordered_set<size_t> mytriples;
//    size_t chkcounter = 0;
//    // __DEBUG

    string chrA;
    string chrB;
    string chrC;
    size_t bpposA = 0;
    size_t bpposB = 0;
    size_t bpposC = 0;

    // ATTENTION! No OMP parallelization here without proper consideration of the exclude set!!!
    for (size_t a = setA.first; a < setA.second; a++) {
        if (excluderange > 0) { // if an exclude range is defined
            chrA = db.getSNPInfo(a).chromosome;
            bpposA = db.getSNPInfo(a).pos_bp;
        }
        // ATTENTION! No OMP parallelization here without proper consideration of the exclude set!!!
        for (size_t b = setB.first; b < setB.second; b++) {
            if (b >= setA.first && a < setB.second && b <= a) { // TODO not good behaviour...
                b = a;
                continue;
            }
            if (excluderange > 0) { // if an exclude range is defined
                chrB = db.getSNPInfo(b).chromosome;
                bpposB = db.getSNPInfo(b).pos_bp;
            }
            atomic<size_t> skipcnt = 0; // TODO we definitely need to calculate skipcnt in advance due to the excluderange!!!
            // apply excluderange between A and B
            if (excluderange == 0 || chrA.compare(chrB) || bpposB - bpposA >= excluderange) {
// TODO parallelization over c makes only sense, if c is large enough!
                omp_set_num_threads(numthreads);
#ifndef DEBUG_SINGLETHREAD
#pragma omp parallel for schedule(dynamic)
#endif
                for (size_t c = setC.first; c < setC.second; c++) {
                    // do not recalculate triples that have already been covered (when permuting a,b,c)
                    if (c >= setA.first && a < setC.second && c <= a) { // TODO not good behaviour...
                        //c = a; does not work for parallel threads!
                        skipcnt++;
                        continue;
                    }
                    if (c >= setB.first && b < setC.second && c <= b) { // TODO not good behaviour...
                        //c = b; does not work for parallel threads!
                        skipcnt++;
                        continue;
                    }
                    if (excluderange > 0) { // if an exclude range is defined
                        chrC = db.getSNPInfo(c).chromosome;
                        bpposC = db.getSNPInfo(c).pos_bp;
                    }
                    if (excluderange == 0 ||
                            ( (chrB.compare(chrC) != 0 || bpposC - bpposB >= excluderange)
                           && (chrA.compare(chrC) != 0 || bpposC - bpposA >= excluderange)
                           // && (chrA.compare(chrB) != 0 || bpposB - bpposA >= excluderange) // this was already checked before snpC loop
                            ) ) {
                        processTriple(a, b, c);
                    }
    //                // DEBUG
    //                mx.lock();
    //                if (a != b && a != c && b != c && mytriples.emplace(xyz(a, b, c)).second)
    //                    chkcounter++;
    ////                else
    ////                    cerr << "XX" << endl;
    //                mx.unlock();
    //                // __DEBUG
                }
            }

            // progress
            numtests_processed += setC.second - setC.first - skipcnt;
            progress.updateProgress(numtests_processed);

        }
    }

//    // DEBUG
//    if (allcounter != chkcounter || allcounter != counter) {
//        cout << "Failed! A:" << setA_.first << "-" << setA_.second << " B:" << setB_.first << "-" << setB_.second << " C:" << setC_.first << "-" << setC_.second << endl;
//        cout << "allcounter:\t" << allcounter << endl;
//        cout << "chkcounter:\t" << chkcounter << endl;
//        cout << "mycounter:\t" << counter << endl;
//    }
//    // __DEBUG

}
