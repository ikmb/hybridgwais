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

#ifndef RESULTHANDLER_H
#define RESULTHANDLER_H

#include <atomic>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <fstream>

#include "hybridsys/Buffer.h"
#include "hybridsys/ThreadUtils.h"

#include "MinMaxHeap.h"
#include "SNPDB.h"
#include "GPUEngine.h"
#include "FPGAHandler.h"
#include "ResultView.h"
#include "ThreadPool.h"
#include "Method.h"
#include "Args.h"
#include "StatusFile.h"

template<class Id = ResultDefaults::default_id_type, class Score = ResultDefaults::default_score_type>
class ResultHandler
{
public:
    explicit ResultHandler(const SNPDB &db, const string &output_file_, int num_threads_, Score threshold_, unsigned long bestResults_, ResultView<Id,Score> &v, const vector<Method> &methods_, Args &args)
        : snpdb(db),
		  output_file(output_file_),
          view(v),
          heaps(num_threads_),
          num_threads(num_threads_),
          threshold(threshold_),
		  ldfilter(args.ldfilter),
          bestResults(bestResults_),
          methods(methods_),
          order(methods_[0].getOrder()),
          debug(args.debug)
    {
        useThreshold = !isnan(threshold);
        useLDFilter = ldfilter > 0;
        useResults = bestResults > 0;

        heaps.resize(num_threads);

		// the last method should be the ld method (if not, i.e. ldfilter is not used, the index points to the last method, which is ok as it will be unused)
        unsigned numscorefields = Method::getScoreFields(methods);
		ldfilter_scoreindex = numscorefields - methods.back().getScoreFields();

		cout << "Results:" << endl;
		if (useThreshold)
			cout << "  Score threshold: " << threshold << endl;
		else
			cout << "  Keeping " << bestResults << " best results." << endl;
		if (useLDFilter)
			cout << "  Applying LD filter threshold: " << ldfilter << endl;
    }

    bool checkID(Result<Id,Score> &r, const FPGAHandler::Block *block, bool apply_correction) {

		// SNPDB local to global ID mapping
		const vector<size_t> &map2global = snpdb.getMappingLtoG();
		size_t id0 = r.getID(0);
		size_t id1 = r.getID(1);
		size_t id2 = order == 3 ? r.getID(2) : 0;

		if (block) { // only used for FPGA processing: FPGA results contain superfluous calculated tests that need to be filtered

			if (id0 > block->snpcount || id1 > block->snpcount || id2 > block->snpcount) {
				stringstream s;
				s << "ResultHandler: Got unallowed (FPGA local) ID: " << id0 << "," << id1;
				if (order == 3)
					s << "," << id2;
				StatusFile::addError(s.str());
				exit(EXIT_FAILURE);
			}

			// the first ID may exceed the required initialization range due to a full PE chain initialization
			if (id0 > block->init_last)
				return false;

			// the second ID may not be in the streaming range due to results generated while PE chain initialization
			if (id1 < block->stream_start)
				return false;

			// TODO for the 3rd ID there's no plan yet...

			const vector<size_t> &idmap = block->idmap;

            // apply excluderange if set
            if (block->excluderange) {
            	if (order == 2) {
					size_t gid0 = map2global[idmap[id0]];
					size_t gid1 = map2global[idmap[id1]];
					string chr0 = snpdb.getSNPInfo(gid0).chromosome;
					int64_t bppos0 = snpdb.getSNPInfo(gid0).pos_bp;
					string chr1 = snpdb.getSNPInfo(gid1).chromosome;
					int64_t bppos1 = snpdb.getSNPInfo(gid1).pos_bp;
					if (chr0.compare(chr1) == 0 && (size_t)abs(bppos1 - bppos0) < block->excluderange)
						return false;
            	} else {
            		// TODO not implemented for third order yet...
            	}
            }

            // swap if necessary
            if (apply_correction && block->swapped) {
            	// TODO no third order yet

            	// swap if outside the overlap region
            	if (id0 < block->stream_start || id1 > block->init_last)
            		swap(id0, id1);
            }

            // map the IDs:
			// FPGA local to SNPDB local ID mapping
			id0 = idmap[id0];
			id1 = idmap[id1];
			id2 = order == 3 ? idmap[id2] : 0;

        } else { // no block provided -> CPU or GPU run -> no need to filter any results
        }

		if (apply_correction) {
			// map to global ID and replace
			r.replaceID(0,map2global[id0]);
			r.replaceID(1,map2global[id1]);
			if (order == 3)
				r.replaceID(2,map2global[id2]);
		}
		return true;
    }

    void processResultView(ResultView<Id,Score> &v, int threadIndex, const FPGAHandler::Block *block = nullptr) {

		// scan through all results
        if (v.getScoreCount()) { // ensure that this method has at least one result
            for(auto it = v.begin(); it != v.end(); ++it) {

                // don't need results if we're working on a permutation
                if(!useResults)
                    continue;

                // filter results by threshold and/or LD
				bool pass = !useThreshold || CmpLess::le(threshold, it.getScore(0)); // if sorting out by threshold is enabled -> compare, else pass.
				pass = pass && (!useLDFilter || // no LD filter
						(CmpLess::le(it.getScore(ldfilter_scoreindex), ldfilter) && (order == 2 || // pairwise LD filter
								(CmpLess::le(it.getScore(ldfilter_scoreindex+1), ldfilter) && CmpLess::le(it.getScore(ldfilter_scoreindex+2), ldfilter))))); // 3rd-order LD filter
				if (pass) {
					// push into heap, maintaining its capacity
					if(heaps[threadIndex].getSize() >= bestResults) {
						if(CmpLess::lt(heaps[threadIndex].getMin(), *it)) {
							Result<Id,Score> r = *it;
							if (checkID(r, block, true)) { // need to check if the ID fits the current sets and correct the result ID immediately
								heaps[threadIndex].push(r);
								heaps[threadIndex].popMin();
							}
						}
					} else {
						Result<Id,Score> r = *it;
						if (checkID(r, block, true)) // need to check if the ID fits the current sets and correct the result ID immediately
							heaps[threadIndex].push(r);
					}
				}
            }
        }
    }

    template<typename Tin, typename Tout>
    void process(typename ThreadPool<Tin, Tout>::inqueue_type &inqueue,
                 typename ThreadPool<Tin, Tout>::outqueue_type &outqueue __attribute__ ((unused)),
                 atomic_bool& termination_request,
                 int threadIndex
                 ) {

        ResultView<Id,Score> v(view); // create a local copy of the view

        ThreadUtils::setThreadName("ResultH:" + to_string(threadIndex));


        while(!termination_request) {
            Tin srcbuf;
            try {
                inqueue.pop(srcbuf);
            } catch(tbb::user_abort &e) { // is thrown whenever the process is cancelled or regularly finished
                break;
            }

            // extract metadata for ID mapping
            const shared_ptr<FPGAHandler::Block> &block = (const shared_ptr<FPGAHandler::Block>&)(srcbuf->getMetadata());

			// set view to current buffer
            v.setBuffer(srcbuf->getData(), (srcbuf->getContentLength()) * v.getResultSize());

            processResultView(v, threadIndex, block.get());
        }
    }

    void flush(const string &cmdline) {

        if(useResults)
            flushResults(cmdline);

    }

    string getResultFileName() const {
        return output_file + ".scores";
    }


private:
    class CmpLess {
    public:
        // returns true if r1.score < r2.score or if the scores are equal, if r1.id < r2.id (comparison started with the first ID component)
        bool operator()(const Result<Id,Score> &r1, const Result<Id,Score> &r2) const {
            return lt(r1,r2);
        }

        static bool lt(const Result<Id,Score> &r1_, const Result<Id,Score> &r2_) {
            const Score &r1 = r1_.getScore(0);
            const Score &r2 = r2_.getScore(0);
            if (r1 < r2 || isnan(r1)) {
                return true;
            } else if (r1 <= r2 && r1 >= r2) { // equal scores
                vector<Id> ids1 = r1_.getID();
                vector<Id> ids2 = r2_.getID();
                for (size_t i = 0; i < ids1.size(); i++) {
                    if (ids1[i] < ids2[i]) {
                        return true;
                    }
                }
            }
            return false;
        }

        static bool le(const Score &r1, const Score &r2) {
            if (r1 <= r2 || isnan(r1))
                return true;
            return false;
        }
    };

//    class CmpEq {
//    public:
//        // returns true if IDs are equal
//        bool operator()(const Result<Id,Score> &r1, const Result<Id,Score> &r2) const {
//            vector<Id> ids1 = r1.getID();
//            vector<Id> ids2 = r2.getID();
//            for (size_t i = 0; i < ids1.size(); i++) {
//                if (ids1[i] != ids2[i])
//                    return false;
//            }
//            return true;
//        }
//    };

    const SNPDB &snpdb;
    string output_file;
    ResultView<Id,Score> view;
//    vector<::MinMaxHeap<Result<Id,Score>,CmpLess,CmpEq>> heaps;
    vector<::MinMaxHeap<Result<Id,Score>,CmpLess>> heaps;

    int num_threads;
    Score threshold;
    Score ldfilter;
    unsigned ldfilter_scoreindex;
    unsigned long bestResults;

    const vector<Method> &methods;
    unsigned order;

    bool useThreshold;
    bool useLDFilter;
    bool useResults;

    bool debug;

    void flushResults(const string &cmdline) {
        // build filename
        stringstream ss;
        ss << output_file << ".scores";
        ofstream results {ss.str()};

        // how many results did we collect?
        unsigned long collected_results = 0ul;
        if (debug)
            cout << "\nResults from " << heaps.size() << " heaps: ";
        for (const ::MinMaxHeap<Result<Id,Score>,CmpLess> &h : heaps) {
        	collected_results += h.getSize();
        	if (debug)
        	    cout << h.getSize() << " ";
        }

        if (debug)
            cout << "\nTotal from heaps: " << collected_results << endl;

        // print command line
        results << "#" << cmdline << endl;

        // print header
        Method::printHeader(methods, results);

        // collect <actual_results> top elements from all our heaps
        unsigned long actual_results = 0;
        while(actual_results < bestResults) {
            auto it = max_element(begin(heaps), end(heaps), [](const ::MinMaxHeap<Result<Id,Score>,CmpLess> &lhs, const ::MinMaxHeap<Result<Id,Score>,CmpLess> &rhs) {
                if(lhs.isEmpty()) return true;
                if(rhs.isEmpty()) return false;
                return CmpLess::lt(lhs.getMax(), rhs.getMax());
            });
            if (it->isEmpty()) // no more elements
                break;

            auto &r = it->getMax();

            // print SNP positions
			for(unsigned id = 0; id < r.getID().size(); id++) {
				const SNPDB::SNPInfo &info = snpdb.getSNPInfo(r.getID(id));
				results << info.chromosome << ":" << info.pos_bp << "\t";
			}
            // print SNP names
            for(unsigned id = 0; id < r.getID().size(); id++) {
				results << snpdb.getSNPInfo(r.getID(id)).variant_id << "\t";
			}
//            // DEBUG
//            if (r.getID(0) == 210458 && r.getID(1) == 566505)
//                cout << "SNP pair has been written to file." << endl;
            // print rest
            results << it->popMax();
            actual_results++;
        }

        results.close();

        StatusFile::addInfo("<p class='pinfo'>Saved " + to_string(actual_results) + " results.</p>");
    }

};

#endif // RESULTHANDLER_H
