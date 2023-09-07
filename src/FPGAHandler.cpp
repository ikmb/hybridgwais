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

#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <tbb/concurrent_queue.h>

#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"

#include "FPGAHandler.h"
#include "FPGAConfigurationGWAIS.h"
#include "Method.h"
#include "utils.h"
#include "StatusFile.h"

using namespace std;

FPGAHandler::FPGAHandler(
        hybridsys::Hybridsys &hybridsys_,
        chrono::milliseconds timeout_,
        hybridsys::BufferFactory<hybridsys::FPGABuffer> &factory,
        const SNPDB &db_,
        const vector<Method> &methods_,
        Progress &progress_,
        const bool &term_request_,
        bool debug_)
    : Processor(db_, methods_, 1, progress_, debug_),
      hybridsys(hybridsys_), timeout(timeout_),
      order(methods_[0].getOrder()), bufferFactory(factory),
      term_request(term_request_),
      debug(debug_)
{

	num_fpgas = hybridsys.getFPGAs().size();

	// skip if there are no FPGAs to use
    if (num_fpgas) {


        fpgaconf.parse(hybridsys.getFPGA(0).getConfigurationData());

        tablesPerBuffer = bufferFactory.getBufferSize() / fpgaconf.getTableSize();

        // initialize DMA channels
        for(hybridsys::FPGA &f : hybridsys.getFPGAs()) { // WARNING! With several FPGAs this could lead to a deadlock if other processes lock FPGAs in another order!
            f.declareDMAChannel(0, hybridsys::FPGA::DMAChannelType::AxiStreamMaster);
            if (debug)
                cout << "FPGA: " << f.getIndex() << " Channel 0: Master" << endl;
            for(unsigned i = 0; i < fpgaconf.getNumChains(); i++) {
                f.declareDMAChannel(i+1, hybridsys::FPGA::DMAChannelType::AxiStreamSlave);
                if (debug)
                    cout << "FPGA: " << f.getIndex() << " Channel " << (i+1) << ": Slave" << endl;
            }
            cout << "Waiting for lock on FPGA #" << f.getIndex() << " (sn: " << f.getSerialNumber() << ")... " << flush;
            f.lock();
            cout << "got it." << endl;
        }

        // prepare mutexes and condition variables for fetchTables threads
        // -> one mtx and cv per thread
        mtxs.resize(num_fpgas);
        cvs.resize(num_fpgas);
        for (size_t f = 0; f < hybridsys.getFPGAs().size(); f++) {
        	for (size_t c = 0; c < fpgaconf.getNumChains(); c++) {
				mtxs[f].push_back(new mutex());
				cvs[f].push_back(new condition_variable());
        	}
        }

        setready_flags.resize(num_fpgas, false);
        buffersRemaining.resize(num_fpgas);
        for (auto &br : buffersRemaining)
        	br.resize(fpgaconf.getNumChains(), 0);
        currentBlocks.resize(num_fpgas, nullptr);
	}

//    processedIntervals.clear();
}

// starts the process to create the contingency tables on the FPGA.
// (loop over initDistributionScheme(), initializeFPGA(), distributeSNPData() for each set combination)
void FPGAHandler::createTables() {
    finished_flag = false;

    // call parent function
    process();

    // The process() method from the parent will subsequently call the processSetPair/Triple() functions below.

    // when we return here, we are finished
    finished_flag = true;

    // to release the fetchTables threads, we need to claim that a new set is ready, but do not touch buffersRemaining[]
    for (unsigned fpgaidx = 0; fpgaidx < num_fpgas; fpgaidx++)
    	setready_flags[fpgaidx] = true;
    for (auto &cv : cvs)
    	for (auto c : cv)
    		c->notify_all();
}


void FPGAHandler::processSetPair(const snprange &setA, const snprange &setB, size_t excluderange) {
	if (term_request)
		return;
	if (SNPDB::isEmpty(setA) || SNPDB::isEmpty(setB))
		return;

	// initialize how to distribute the current pair of sets
	vector<shared_ptr<Block>> distribution = initDistribution(setA, setB, excluderange);

	// build a concurrent queue from the distribution and run one thread per available FPGA
	tbb::concurrent_queue<shared_ptr<Block>> distrqueue;
	for (const auto &d : distribution)
		distrqueue.push(d);

	// parallelize initialization over the number of available FPGAs:
	// the distribution was created to have num_fpgas blocks of nearly equal size next to each other
	omp_set_num_threads(num_fpgas);
#pragma omp parallel for
	for (unsigned fpgaidx = 0; fpgaidx < num_fpgas; fpgaidx++) {
		shared_ptr<Block> blockptr;
		while (distrqueue.try_pop(blockptr)) { // all items are already in the queue. if the queue is empty, the thread can terminate

			// acquire lock on all mutexes before initializing FPGA with this pair of sets
			// -> this ensures all data from previous set is fetched
			// Note that each thread uses its own mux such that we should not run into a deadlock situation here.
			for (auto m : mtxs[fpgaidx])
				m->lock();

			// set current block for current FPGA
			currentBlocks[fpgaidx] = blockptr;
			buffersRemaining[fpgaidx] = blockptr->buffercount;

			initializeFPGA(fpgaidx);
			setready_flags[fpgaidx] = true;

			// the fetchTables threads are waiting each on their cv, we can start them now!
			// -> we notify all threads and unlock our mutexes
			for (auto cv : cvs[fpgaidx])
				cv->notify_all();

			for (auto m : mtxs[fpgaidx])
				m->unlock();

			// before we reset the setready_flag, we need to wait until all results are fetched
			while (!isBlockFinished(fpgaidx) && !term_request)
				this_thread::yield();

			setready_flags[fpgaidx] = false;
		}
	}
}

void FPGAHandler::processSetTriple(const snprange &setA __attribute__((unused)), const snprange &setB __attribute__((unused)), const snprange &setC __attribute__((unused)), size_t excluderange_ __attribute__((unused))) {
    if (term_request)
        return;
    if (SNPDB::isEmpty(setA) || SNPDB::isEmpty(setB) || SNPDB::isEmpty(setC))
		return;
    // 3rd order FPGA support is deprecated
}

vector<shared_ptr<FPGAHandler::Block>> FPGAHandler::initDistribution(const snprange &setA_, const snprange &setB_, size_t excluderange) {
	// This method is intended for 2nd order tests.
	if (order > 2)
		return vector<shared_ptr<FPGAHandler::Block>>();

	// Note: both provided sets are not empty!

	// fill a single block with no restrictions first
	shared_ptr<Block> initblockptr(new Block());
	{
    	Block &block = *initblockptr;

    	vector<snprange> &sets = block.sets;

    	// It is necessary that the provided sets in the block
    	// DO NOT OVERLAP!!!

    	// for more efficiency, the smaller set should be the first set
    	bool swapped = SNPDB::lengthOf(setA_) > SNPDB::lengthOf(setB_); // need to swap?
    	snprange setA = swapped ? setB_ : setA_;
    	snprange setB = swapped ? setA_ : setB_;
    	block.swapped = swapped;

    	// we split the provided sets in three intervals:
    	// 1. the overlapping region
    	// 2a. region of A before overlap
    	// 2b. region of B before overlap //( -> will be empty after swapping anyway if there is an overlap )
    	// 3a. region of A after overlap
    	// 3b. region of B after overlap
    	// Note that either 2a or 2b is empty, same for 3a and 3b.
    	snprange ov = SNPDB::intersec(setA, setB); // 1.
    	pair<snprange,snprange> cuta = SNPDB::cut(setA, ov); // 2a.,3a.
    	pair<snprange,snprange> cutb = SNPDB::cut(setB, ov); // 2b.,3b.

    	if ((SNPDB::isEmpty(cuta.first) && SNPDB::isEmpty(cuta.second)) || (SNPDB::isEmpty(cutb.first) && SNPDB::isEmpty(cutb.second))) { // one set is completely included in the other
			// special treatment where we first push the overlapping part followed by the non-overlapping parts
			sets.push_back(ov);
			if (!SNPDB::isEmpty(cuta.first))
				sets.push_back(cuta.first);
			else if (!SNPDB::isEmpty(cutb.first))
				sets.push_back(cutb.first);
			if (!SNPDB::isEmpty(cuta.second))
				sets.push_back(cuta.second);
			else if (!SNPDB::isEmpty(cutb.second))
				sets.push_back(cutb.second);

			block.init_last = SNPDB::lengthOf(ov) - 1; // note that ov is not empty and we need FPGA local IDs here
			block.stream_start = 0; // all pairs until init_last are required

			// need to swap if the streaming part belongs to A
			if (!SNPDB::isEmpty(cuta.first) || !SNPDB::isEmpty(cuta.second))
				block.swapped = !block.swapped;
		} else {
			snprange apart;
			if (!SNPDB::isEmpty(cuta.first))
				apart = cuta.first;
			else // this means cuta.second is not empty, otherwise it would be handled by the case above
				apart = cuta.second;
			sets.push_back(apart);

			if (!SNPDB::isEmpty(ov)) // otherwise sets A and B would be distinct -> no problem here!
				sets.push_back(ov);

			if (!SNPDB::isEmpty(cutb.first))
				sets.push_back(cutb.first);
			else // this means cutb.second is not empty, otherwise it would be handled by the case above
				sets.push_back(cutb.second);

			block.init_last = SNPDB::lengthOf(apart) + SNPDB::lengthOf(ov) - 1; // note that we need FPGA local IDs here
			block.stream_start = SNPDB::lengthOf(apart); // cannot be empty, stream over ov and cutb
		}

        block.snpcount = 0;
        for (const snprange &set : sets)
            block.snpcount += SNPDB::lengthOf(set);
	}


	// distribute the initial block over the number of FPGAs with an approximately equal amount of tests for each FPGA (same area)
	vector<shared_ptr<Block>> fpgadistribution;
	// inefficient FPGA usage if the number of init SNPs is below the length of the PE chain
	size_t chainlength = fpgaconf.getNumPEPerChain() * fpgaconf.getNumChains();
	if (num_fpgas > 1 && initblockptr->init_last > chainlength) {
		// NOTE: we are expecting only a small number of FPGAs here, so we don't check if the parts get too small
		const Block &initblock = *initblockptr;

		// calculate area to distribute
		size_t area = (initblock.snpcount - initblock.init_last - 1) * (initblock.init_last + 1) // rectangle right of init triangle
				    + ((initblock.init_last+1)*initblock.init_last) / 2 // init_triangle
					- (initblock.stream_start*(initblock.stream_start-1)) / 2; // minus triangle before stream_start
//		// DEBUG
//		cout << "total area: " << area << endl;

		area = area / num_fpgas;

//		// DEBUG
//		cout << "area per FPGA: " << area << endl;

		// calculate for all FPGAs the number of SNPs for the init part (called "width")
		vector<size_t> w(num_fpgas, 0);
		{
			size_t remainder = initblock.init_last+1; // number of remaining SNPs of init block to distribute
			size_t snpcount = initblock.snpcount;
			size_t stream_start = initblock.stream_start;
			for (unsigned f = 0; f < num_fpgas-1; f++) { // don't calculate for the last FPGA -> simply take the remainder
				// 2 cases: either w[f] will be less than stream_start, then the area is a simple rectangle,
				// or otherwise, we need to cut a triangle from the upper left corner
				w[f] = area / (snpcount - stream_start); // rectangle case
				if (w[f] > stream_start) { // need to calculate for case 2 (w[f] == stream_start will still be a rectangle!)
					// p,q equation, second solution with subtracting the square root, as width may not exceed snpcount
					w[f] = (size_t)round((double)snpcount - 0.5 - sqrt(snpcount*(snpcount-1) - stream_start*(stream_start-1) - 2.0*area - 0.25));
				}
				// need to adjust the remainder and the remaining block now
				remainder -= w[f];
				snpcount -= w[f];
				stream_start = w[f] > stream_start ? 0 : (stream_start - w[f]);
			}
			// set the remainder as width for the last FPGA
			w.back() = remainder;

//			// DEBUG
//			cout << "final widths: ";
//			for (auto x:w)
//				cout << " " << x;
//			cout << endl;
		}

		// depending on the calculated width, create new blocks for the distribution
		{
			size_t remainder = initblock.init_last+1; // number of remaining SNPs to distribute
			size_t snpcount = initblock.snpcount;
			size_t stream_start = initblock.stream_start;
			vector<snprange> remaining_sets = initblock.sets;
			// formally, we also need to check if the calculated width is greater than zero
			for (unsigned f = 0; f < num_fpgas && w[f] > 0; f++) {
				fpgadistribution.emplace_back(new Block());
				Block &block = *fpgadistribution.back();
				block.swapped = initblock.swapped;
				if (w[f] < stream_start) { // simple rectangle
					block.init_last = w[f]-1;
					block.stream_start = w[f];
					// need to remove the SNPs from new stream_start to old stream_start
					block.snpcount = snpcount - (stream_start - w[f]);
					block.sets = SNPDB::extract(remaining_sets, 0, w[f]);
					vector<snprange> moresets = SNPDB::extract(remaining_sets, stream_start, snpcount-stream_start);
					block.sets.insert(block.sets.end(), moresets.begin(), moresets.end());
				} else { // triangle + rectangle
					block.init_last = w[f]-1;
					block.stream_start = stream_start; // if positive, the first edge of the triangle is cut off
					block.snpcount = snpcount;
					block.sets = remaining_sets;
				}
				// set remainders
				remaining_sets = SNPDB::extract(remaining_sets, w[f], snpcount-w[f]);
				remainder -= w[f];
				snpcount -= w[f];
				stream_start = w[f] < stream_start ? (stream_start - w[f]) : 0;
			}
		}

//		// DEBUG
//		cout << "distributed areas: " << endl;
//		for (const auto &d : fpgadistribution) {
//			size_t t = (d->init_last+1)*(d->init_last)/2;
//			size_t tm = (d->stream_start)*(d->stream_start-1)/2;
//			size_t r = (d->snpcount - d->init_last - 1) * (d->init_last+1);
//			cout << "  t = " << (d->init_last+1) << " x " << (d->init_last) << " / 2 = " << t << endl;
//			if (d->stream_start)
//				cout << "  tm = " << (d->stream_start) << " x " << (d->stream_start-1) << " / 2 = " << tm << endl;
//			cout << "  r = " << (d->snpcount - d->init_last - 1) << " x " << (d->init_last+1) << " = " << r << endl;
//			if (d->stream_start)
//				cout << "  t - tm + r = " << (t - tm + r) << endl;
//			else
//				cout << "  t + r = " << (t+r) << endl;
//		}
	} else // one FPGA or init block too small
		fpgadistribution.push_back(initblockptr);


	// check if the size of each block in the current FPGA distribution meets the requirements and create the final distribution
	vector<shared_ptr<Block>> distribution;
	for (size_t fidx = 0; fidx < fpgadistribution.size(); fidx++) {
		shared_ptr<Block> blockptr = fpgadistribution[fidx]; // we make a copy here as we might clean the distribution below (in the case we need to split up due to memory issues)
		Block &block = *blockptr;
		size_t samplecnt = db.getSampleCount();
		size_t maxsnps = fpgaconf.getMaxGenotypeCount() / samplecnt;
		if (block.snpcount > maxsnps) { // too large! need to distribute!

			size_t ovsize = block.init_last+1;
			size_t novsize = block.snpcount - ovsize;
			// number of parts necessary for the ov part: Two(!) distinct parts need to fit into memory
			size_t num_ovp = divideRounded(ovsize, maxsnps/2);
			size_t ovpsize = ovsize / num_ovp; // note that the remainder ovsize % num_ovp must not be forgotten!
			// the size of a non-ov part (after init_last) must not exceed the size of an ov part, implicating the number of non-ov parts
			size_t num_novp = max(1l, divideRounded(novsize, ovpsize)); // corner case: we need at least one block, even if it's of size 0.
			size_t novpsize = novsize / num_novp; // due to rounding, this might be slightly different to ovpsize. do not forget the remainder as well!

//			// DEBUG
//			cout << "  ovsize:   " << ovsize << endl;
//			cout << "  novsize:  " << novsize << endl;
//			cout << "  num_ovp:  " << num_ovp << endl;
//			cout << "  num_novp: " << num_novp << endl;
//			cout << "  ovpsize:  " << ovpsize << endl;
//			cout << "  novpsize:  " << novpsize << endl;

			// generate all combinations of parts now for the distribution
			size_t curr_posa = 0;
			for (size_t ovi = 0; ovi < num_ovp; ovi++) {
				size_t curr_sizea = ovpsize + (ovi < (ovsize % num_ovp) ? 1 : 0); // distribute the remainder!
				size_t curr_posb = curr_posa + curr_sizea;
				for (size_t novi = 0; novi < num_novp + num_ovp-ovi-1; novi++) { // we enter at least once, even if the size of the nov part is zero, and we need to handle the extra "squares" in the triangle before non-ov region
					size_t curr_sizeb;
					if (novi < num_ovp-ovi-1) { // this is an extra "square" in the triangle before non-ov
						curr_sizeb = ovpsize + ((ovi+novi+1) < (ovsize % num_ovp) ? 1 : 0); // distribute the remainder, but for the corresponding square of ovsize!
					} else {
						curr_sizeb = novpsize + ((novi - (num_ovp-ovi-1)) < (novsize % num_novp) ? 1 : 0); // distribute the remainder!
					}

//					// DEBUG
//					cout << ovi << "/" << novi << ": posa: " << curr_posa << " sizea: " << curr_sizea << " posb: " << curr_posb << " sizeb: " << curr_sizeb << endl;

					// only continue if the block will be valid, i.e. we must be the first triangle or have a valid stream range AND the original stream start must be before or in the new stream region
					if ((novi == 0 || curr_sizeb > 0) && block.stream_start < curr_posb + curr_sizeb) {
						// create a new sub block for the current combination
						shared_ptr<Block> newbptr(new Block());
						Block &newblock = *newbptr;
						newblock.init_last = curr_sizea - 1;
						size_t stream_start_correction = 0;
						if (block.stream_start <= curr_posb) { // complete or part of triangle, but complete rectangle
							// NOTE: novi == 0 means curr_posb = curr_posa + curr_sizea, so b begins directly after a with no gap
							newblock.stream_start = novi == 0 ? (block.stream_start > curr_posa ? block.stream_start - curr_posa : 0) : curr_sizea; // consider the triangle only at the beginning of a line, but cut the beginning according to the original stream_start then!
						} else { // no triangle, only part of the rectangle
							newblock.stream_start = curr_sizea; // the part of the rectangle is cut by reducing the number of SNPs, so streaming still begins directly after part a, but independent of novi
							// need to adjust snpcount and the setsB by this number
							stream_start_correction = block.stream_start - curr_posb;
						}
						newblock.snpcount = curr_sizea + curr_sizeb - stream_start_correction;
						newblock.sets = SNPDB::extract(block.sets, curr_posa, curr_sizea); // the required init set(s)
						vector<snprange> setsB = SNPDB::extract(block.sets, curr_posb + stream_start_correction, curr_sizeb - stream_start_correction); // the stream set(s)
						if (!setsB.empty())
							newblock.sets.insert(newblock.sets.end(), setsB.begin(), setsB.end());
						newblock.swapped = novi < num_ovp-ovi-1 ? false : block.swapped; // do not swap in the original ov region, otherwise swap if the original block should be swapped

						// insert
						distribution.push_back(newbptr);
					}
					curr_posb += curr_sizeb;
				} // end for nov parts
				curr_posa += curr_sizea;
			} // end for ov parts
		} else // not too large -> simply take over the current block
			distribution.push_back(blockptr);
	} // end block


	// for all final distribution blocks, fill the remaining fields
	for (shared_ptr<Block> &block : distribution) {
        vector<size_t> &idmap = block->idmap;
        idmap.resize(block->snpcount);
        size_t currid = 0;
        for (const snprange &set : block->sets)
            for (size_t i = set.first; i < set.second; i++)
                idmap[currid++] = i;

        block->excluderange = excluderange;

        // calculate expected result counters
		vector<size_t> resultsPerChannel;
		resultsPerChannel.resize(fpgaconf.getNumChains());
		block->buffercount.resize(fpgaconf.getNumChains());
		block->residualresults.resize(fpgaconf.getNumChains());
		size_t resultcount_filtered = 0;

		unsigned chain;
		size_t num_pes = fpgaconf.getNumPEPerChain() * fpgaconf.getNumChains();
		for(size_t snp = 0; snp < min((size_t)increaseToMultipleInt(block->init_last+1, num_pes), block->snpcount); snp++) { // FPGA computes blocks of num_pes SNPs in each round, so, also SNPs after init_last are computed
			chain = (snp % num_pes) / fpgaconf.getNumPEPerChain(); // for two channels, this results in either 0 or 1 (indicating DMA channels 1 or 2 resp.)
			resultsPerChannel[chain] += (block->snpcount - snp - 1);
			size_t round_end = increaseToMultipleInt(snp+1, num_pes); // exclusive
			if (round_end < block->stream_start) // FPGA jumps over the gap between the last SNP in the round and the stream start
				resultsPerChannel[chain] -= block->stream_start - round_end;
			// calculate number of results after filtering
			if (snp <= block->init_last) { // only SNPs inside initialization area generate real (unfiltered) results
				resultcount_filtered += block->snpcount - max(snp, block->stream_start) - 1; // only SNPs starting from stream start are recognized
			}
		}

		size_t totalbuffers = 0;
		for(unsigned i = 0; i < fpgaconf.getNumChains(); i++) {
			size_t buffer_count = resultsPerChannel[i] / tablesPerBuffer + (resultsPerChannel[i] % tablesPerBuffer ? 1 : 0);
			block->buffercount[i] = buffer_count;
			totalbuffers += buffer_count;
			block->residualresults[i] = resultsPerChannel[i] % tablesPerBuffer;
			if (block->residualresults[i] == 0ull && block->buffercount[i] > 0ull)
				block->residualresults[i] = tablesPerBuffer;
		}

		block->resultcount_filtered_per_buffer = resultcount_filtered / totalbuffers;

    }


	if (debug) {
    	for (size_t b = 0; b < distribution.size(); b++) {
    		const Block &block = *distribution[b];
			cout << "Distribution block " << b << ":" << endl;
			for (const auto &s : block.sets)
				cout << " " << SNPDB::printSNPRange(s) << endl;
			cout << " snpcount:     " << block.snpcount << endl;
			cout << " init_last:    " << block.init_last << endl;
			cout << " stream_start: " << block.stream_start << endl;
			cout << " swapped:      " << (block.swapped ? "yes" : "no") << endl;
			cout << " expected tables per chain:" << endl;
			for(unsigned i = 0; i < fpgaconf.getNumChains(); i++) {
				size_t resultsperchannel = block.buffercount[i] > 0 ? (block.buffercount[i]-1)*tablesPerBuffer + block.residualresults[i] : 0;
				cout << "  DMA channel " << (i+1) << ": " << resultsperchannel << " tables, ";
				cout << block.buffercount[i] << " buffers, ";
				cout << si_binary(block.buffercount[i] * bufferFactory.getBufferSize());
				cout << ", " << block.residualresults[i] << " tables in last buffer." << endl;
			}
			cout << " average filtered result count per buffer: " << block.resultcount_filtered_per_buffer << endl;
    	}
    }

	return distribution;
}


// initialize constants for FPGA run according to current distribution scheme, FPGA needs to be idle when calling this function
void FPGAHandler::initializeFPGA(unsigned fpgaidx) {

	const shared_ptr<Block> &block = currentBlocks[fpgaidx];

    // initialize memory for constants
    auto buffer = bufferFactory.get(true);
    if(buffer->getSize() < num_constants * sizeof(uint64_t)) {
        StatusFile::addError("Buffer is not large enough to hold the FPGA initialization constants (need at least " + to_string(num_constants * sizeof(uint64_t)) + " bytes).");
        exit(EXIT_FAILURE);
    }

    uint64_t *constants = reinterpret_cast<uint64_t*>(buffer->getData());

    // sync word
    constants[0] = 0xDEADBEEFDEADBEEFULL;
    constants[1] = 0xDEADBEEFDEADBEEFULL;
    constants[2] = 0xDEADBEEFDEADBEEFULL;
    constants[3] = 0xDEADBEEFDEADBEEFULL;

    // data set case/control/SNP counters
    size_t scp = db.getSampleCountPadded();
    constants[4] = (db.getCaseCountPadded() << 32) | scp;
    constants[5] = block->snpcount; // local SNP count for FPGA
    // NOTE: for a complete analysis set init last to num_snps-1 and stream start to 0
    constants[6] = (block->stream_start << 32) | block->init_last; // lower 32bit: init last index (inclusive); upper 32bit: stream start idx
    // help the FPGA to calculate the stream start address and the round address offset:
    size_t round_addr_offset = fpgaconf.getNumPEPerChain() * fpgaconf.getNumChains() * 8 * (scp/256 + (scp%256 ? 1 : 0)); // NUM_PE * LONGS_PER_RAMWORD * num_sample_blocks
    size_t stream_start_addr = block->stream_start * 8 * (scp/256 + (scp%256 ? 1 : 0)); // stream_start_idx * LONGS_PER_RAMWORD * num_sample_blocks
    constants[7] = (round_addr_offset << 32) | stream_start_addr;

    // actual data set sizes and buffer information
    constants[8] = ((db.getSNPSize() / 32) * block->snpcount - 1); // last transmission word for SNP data
    constants[9] = (bufferFactory.getBufferSize() / fpgaconf.getTableSize()) << 32;
    constants[9] |= bufferFactory.getBufferSize() / (256ULL/8);
    constants[10] = 0;
    constants[11] = 0;

    // send constants
    hybridsys::FPGA& f = hybridsys.getFPGA(fpgaidx);
	f.writeDMA(*buffer, 0, timeout, num_constants * sizeof(uint64_t));

    // send SNP data
    sendSNPData(fpgaidx);
}


/**
 * @brief Distribute SNP data across FPGAs. After distribution,the ctable creation starts immediately.
 */
void FPGAHandler::sendSNPData(unsigned fpgaidx) {

    auto buffer = bufferFactory.get();
    const shared_ptr<Block> &block = currentBlocks[fpgaidx];

	for (const snprange &set : block->sets) {
		const unsigned char *start = db[set.first];
		size_t length = SNPDB::isEmpty(set) ? 0 : (set.second - set.first) * db.getSNPSize();
		if (debug) {
			cout << "SNP interval & length (bytes) for FPGA " << fpgaidx << ": [" << set.first << "," << set.second << ") " << length << endl;
		}
		while(length > 0) {
			const size_t curr_length = min(length, buffer->getSize());
			memcpy(buffer->getData(), start, curr_length);
			hybridsys.getFPGA(fpgaidx).writeDMA(*buffer, 0, timeout, curr_length);
			length -= curr_length;
			start += curr_length;
		}
	}

}


void FPGAHandler::fetchTables(tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &inqueue __attribute__((unused)),
             tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> &outqueue,
             atomic_bool& termination_request,
             int threadIndex
             ) {

    unsigned fpga_index = threadIndex / fpgaconf.getNumChains();
    unsigned channel_index = threadIndex % fpgaconf.getNumChains(); // 1 or 2 mapped to 0 or 1 resp.

    hybridsys::FPGA &fpga = hybridsys.getFPGAs()[fpga_index];

    ThreadUtils::setThreadName("FPGA:" + to_string(fpga.getIndex()) + "," + to_string(channel_index+1));

    while (!finished_flag) {

        unique_lock<mutex> lock(*(mtxs[fpga_index][channel_index]));
        while ((!setready_flags[fpga_index] || buffersRemaining[fpga_index][channel_index] == 0) && !finished_flag) {
//            // DEBUG
//            cerr << "fetch thr " << threadIndex << " sleeping" << endl;

            cvs[fpga_index][channel_index]->wait(lock);

//            // DEBUG
//            cerr << "woke up " << threadIndex << endl;
        }

        bool done = buffersRemaining[fpga_index][channel_index] == 0;

        while(!termination_request && !done) {

            auto b = bufferFactory.get();

//            // DEBUG
//            cerr << "waiting for buffer " << threadIndex << endl;

            hybridsys::FPGA::Result result = fpga.readDMA(*b, channel_index+1, timeout);
            switch(result) {
            case hybridsys::FPGA::Result::Success:
                if (buffersRemaining[fpga_index][channel_index] == 1) // last buffer
                    b->setContentLength(currentBlocks[fpga_index]->residualresults[channel_index]);
                else
                    b->setContentLength(tablesPerBuffer);

                {
                    // add current distribution block as metadata to the current buffer
                    b->setMetadata((shared_ptr<void>)currentBlocks[fpga_index]);
                }

//                // DEBUG
//                {
//					// dump the buffer
//					ofstream o(string("tmp.buf").append(to_string(channel_index)));
//					uint32_t *buf = reinterpret_cast<uint32_t*>(b->getData());
//					// print the first ctables in the buffer
//	//				for (int j = 0; j < b->getContentLength(); j++) {
//	//				for (int j = 0; j < 2000; j++) {
//					for (int j = 0; j < b->getSize()/44; j++) {
//						o << dec << buf[10] << "\t" << (buf[9]+1) << " : " << hex;
//						for (int i = 0; i < 9; i++) {
//							o << "\t" << setw(8) << setfill('0') << *buf++;
//						}
//						o << "\n";
//						buf += 2; // ID was already printed
//					}
//					o.close();
//                }
//                // __DEBUG

                progress_mx.lock();
                numtests_processed += currentBlocks[fpga_index]->resultcount_filtered_per_buffer; // this reflects the average number of tests per buffer for this FPGA (with a small rounding difference)
                progress.updateProgress(numtests_processed);
                progress_mx.unlock();

                outqueue.push(b);
                buffersRemaining[fpga_index][channel_index]--;

//                // DEBUG
//                cerr << "processed buffer " << threadIndex << endl;

                break;

            case hybridsys::FPGA::Result::Timeout:
                cerr << "Timeout on channel " << (channel_index+1) << "!" << endl;
                timeout_flag = true; // this assignment is not thread-safe, but it's also not critical since it will never be deasserted anyway
                done = true;
                break;

            case hybridsys::FPGA::Result::Cancelled:
                done = true;
                break;
            }
            done |= buffersRemaining[fpga_index][channel_index] == 0;

        } // end while (!done)

    } // end while(!finished_flag)

}
