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

#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cassert>

// for shared memory access:
extern "C" {
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
}

#include "utils.h"
#include "StatusFile.h"

#include "SNPDB.h"

using namespace std;

SNPDB::SNPDB(unsigned minCases, unsigned minControls, bool debug_)
    : order(0),
	  num_snps_global(0),
      num_snps_local(0),
      num_cases(0),
      num_controls(0),
      num_cases_padded(0),
      num_controls_padded(0),
      word_size(1),
      snp_size(0),
      case_size(0),
      current_case_offset {nullptr, 0},
      current_ctrl_offset {nullptr, 0},
      current_init_snp(0),
      allow_buffer_resize(false),
      buffer_size(0),
      buffer(nullptr),
      min_cases(minCases),
      min_controls(minControls),
	  debug(debug_)
{
}

SNPDB::~SNPDB() {
    if(buffer && buffer_fd == -1) {
        delete[] buffer;
    }

    if(buffer && buffer_fd != -1) {
        munmap(buffer, buffer_size);
        close(buffer_fd);
        buffer_fd = -1;
    }

    buffer = nullptr;
}

void SNPDB::initialize(size_t word_size_) {

    this->word_size = word_size_;

    // parse and apply user regions:
    applyRegions();

    // apply padding to actual case/control counts
    num_cases_padded = increaseToMultipleInt(max(num_cases, min_cases), 16ul);
    num_controls_padded = increaseToMultipleInt(max(num_controls, min_controls), 16ul);

    // count bytes required for cases and controls (4gts per byte)
    case_size = num_cases_padded / 4;
    size_t ctrl_size = num_controls_padded / 4;

    snp_size = case_size + ctrl_size;

    // round up to the next word_size
    if(snp_size % word_size)
        snp_size += word_size - (snp_size % word_size);

    buffer_size = snp_size * num_snps_local;

    buffer = new unsigned char[buffer_size];
    memset(buffer, buffer_init_value, buffer_size);
    allow_buffer_resize = false;

    resetCurrOffsets();

}

void SNPDB::applyRegions() {
    // 1. parse region strings
    snprange sA  = parseSetRange(usersetA);
    snprange sB  = parseSetRange(usersetB);
    snprange sC  = parseSetRange(usersetC);
    snprange exs = parseSetRange(userexset);
    region.setA_g  = sA;
    region.setB_g  = sB;
    region.setC_g  = sC;
    region.exset_g = exs;

    // 2. create flags for ignored SNPs
    vector<snprange> sets;
    snprange empty = make_pair(0,0); // will be used later with explicit borders for empty intervals
    auto sAcutex = cut(sA, exs);
    auto sBcutex = region.useexsetA ? make_pair(sB, empty) : cut(sB, exs);
    auto sCcutex = region.useexsetA ? make_pair(sC, empty) : cut(sC, exs);
    sets.push_back(sAcutex.first);
    sets.push_back(sAcutex.second);
    sets.push_back(sBcutex.first);
    sets.push_back(sBcutex.second);
    sets.push_back(sCcutex.first); // empty if not a third order test!
    sets.push_back(sCcutex.second);

//    // DEBUG
//    cout << "\n cutexA: " << printSNPRange(sAcutex.first) << " + " << printSNPRange(sAcutex.second) << endl;
//    cout << " cutexB: " << printSNPRange(sBcutex.first) << " + " << printSNPRange(sBcutex.second) << endl;
//    cout << " cutexC: " << printSNPRange(sCcutex.first) << " + " << printSNPRange(sCcutex.second) << endl;

    // create union of all sets to determine required sets
    auto unit = unite(sets);

//    // DEBUG
//    cout << "\nUnited: " << endl;
//    for (const auto &u : unit)
//    	cout << " " << printSNPRange(u) << endl;
//    // __DEBUG

    // swipe through sets to generate the ignored flags
    ignored_snps.clear();
    ignored_snps.reserve(num_snps_global);
    size_t curr = 0;
    num_snps_local = 0;

    for (const auto &u : unit) {
        num_snps_local += u.second - u.first;
        if (u.first - curr > 0)
            ignored_snps.insert(ignored_snps.end(), u.first - curr, true); // interval before current range is ignored
        ignored_snps.insert(ignored_snps.end(), u.second-u.first, false); // don't ignore specified range
        curr = u.second;
    }
    if (curr < num_snps_global) // ignore remainder after last range
        ignored_snps.insert(ignored_snps.end(), num_snps_global-curr, true);

    // 3. generate map for local SNPs to global index
    mapLtoG.clear();
    mapLtoG.reserve(num_snps_local);
    for (const auto &u : unit) {
        for (size_t l = u.first; l < u.second; l++)
            mapLtoG.push_back(l);
    }

    // 4. generate local ranges

    // map ranges above to local ranges
    unordered_map<size_t,size_t> mapGtoL;
    sort(sets); // sets are now ordered and empties removed

    if (debug) {
		cout << "\nSorted sets:" << endl;
		for (const auto &s : sets)
			cout << " " << printSNPRange(s) << endl;
    }

    snprange highestrange = make_pair(0,0);
    for (const auto &s : sets) {
        if (isEmpty(highestrange)) { // must be the first set
            mapGtoL[s.first] = 0;
        } else if (s.first < highestrange.second) { // overlap
            // add distance to beginning of first interval to the mapping of the last range
            mapGtoL[s.first] = s.first - highestrange.first + mapGtoL[highestrange.first];
        } else { // gap (distance >= 0) between this and last set -> is removed for local mapping
            mapGtoL[s.first] = mapGtoL[highestrange.second]; // in the case where the distance is 0, the values are equal, but this is ok
        }
        // this is common for all new entered sets
        mapGtoL[s.second] = s.second - s.first + mapGtoL[s.first];
        if (s.second > highestrange.second)
        	highestrange = s;
    }
    pair<snprange,snprange> sAmapped;
    pair<snprange,snprange> sBmapped;
    pair<snprange,snprange> sCmapped;
    sAmapped.first = make_pair(mapGtoL[sAcutex.first.first], mapGtoL[sAcutex.first.second]);
    if (!isEmpty(sAcutex.second))
        sAmapped.second = make_pair(mapGtoL[sAcutex.second.first], mapGtoL[sAcutex.second.second]);
    sBmapped.first = make_pair(mapGtoL[sBcutex.first.first], mapGtoL[sBcutex.first.second]);
    if (!isEmpty(sBcutex.second))
        sBmapped.second = make_pair(mapGtoL[sBcutex.second.first], mapGtoL[sBcutex.second.second]);
    if (!isEmpty(sCcutex.first))
        sCmapped.first = make_pair(mapGtoL[sCcutex.first.first], mapGtoL[sCcutex.first.second]);
    if (!isEmpty(sCcutex.second))
        sCmapped.second = make_pair(mapGtoL[sCcutex.second.first], mapGtoL[sCcutex.second.second]);

//    //DEBUG
//    cout << "Mapped:" << endl;
//    cout << " " << printSNPRange(sAmapped.first) << " + " << printSNPRange(sAmapped.second) << endl;
//    cout << " " << printSNPRange(sBmapped.first) << " + " << printSNPRange(sBmapped.second) << endl;
//    cout << " " << printSNPRange(sCmapped.first) << " + " << printSNPRange(sCmapped.second) << endl;

    // combine local intervals
    pair<snprange,snprange> sAlocal;
    snprange sBlocal;
    snprange sClocal;
    if (!isEmpty(sAmapped.second) && sAmapped.first.second == sAmapped.second.first) { // equal start and end -> combine
        sAlocal.first = make_pair(sAmapped.first.first, sAmapped.second.second);
    } else // nothing to do -> simply copy
        sAlocal = sAmapped;
    // per definition, if there exists a non-empty second in B or C, they must be mergeable (as the user exclude set was applied equally to all sets then)
    if (!isEmpty(sBmapped.second)) {
        assert(sBmapped.first.second == sBmapped.second.first); // there must be equal start and end here! -> combine
        sBlocal = make_pair(sBmapped.first.first, sBmapped.second.second);
    } else // simply copy
        sBlocal = sBmapped.first;
    if (!isEmpty(sCmapped.second)) {
        assert(sCmapped.first.second == sCmapped.second.first); // there must be equal start and end here! -> combine
        sClocal = make_pair(sCmapped.first.first, sCmapped.second.second);
    } else // simply copy
        sClocal = sCmapped.first; // could also be empty here

//    //DEBUG
//	cout << "Local mapped:" << endl;
//	cout << " " << printSNPRange(sAlocal.first) << " + " << printSNPRange(sAlocal.second) << endl;
//	cout << " " << printSNPRange(sBlocal) << endl;
//	cout << " " << printSNPRange(sClocal) << endl;

    empty = make_pair(0,0); // reset "empty" interval with consistent borders

    // determine first interval (part left of first overlap of A and B):
    // in the case of 3way, if the exclude range is only for A and not completely covered by B,
    // we need special treatment
    if (isLeftOf(sAlocal.first, sBlocal)) { // A starts left of B
        if (sAlocal.first.second < sBlocal.first) { // gap between A and B
            region.l0 = sAlocal.first;
            // gap is completely left of B?
            if (!isEmpty(sAlocal.second) && isLeftOf(sAlocal.second, sBlocal)) {
                if (sAlocal.second.second < sBlocal.first) // A finishes before B
                    region.l0b = sAlocal.second; // is not empty!
                else
                    region.l0b = make_pair(sAlocal.second.first, sBlocal.first); // is not empty!
            }
        } else // no gap between A and B
            region.l0 = make_pair(sAlocal.first.first, sBlocal.first);
        region.l0isA = true;
    } else if (isLeftOf(sBlocal, sAlocal.first)) { // A starts right of B
        if (sBlocal.second < sAlocal.first.first) // gap between A and B
            region.l0 = sBlocal;
        else // no gap between A and B
            region.l0 = make_pair(sBlocal.first, sAlocal.first.first);
        region.l0isA = false;
    } else
        region.l0 = empty;
    if (SNPDB::isEmpty(region.l0b))
    	region.l0b = make_pair(region.l0.second, region.l0.second); // keep interval borders consistent
    empty = make_pair(region.l0b.second, region.l0b.second); // keep consistent interval borders

    // determine first overlap of A and B
    region.l1 = intersec(sAlocal.first, sBlocal);
    if (SNPDB::isEmpty(region.l1))
    	region.l1 = empty; // keep consistent borders!
    empty = make_pair(region.l1.second, region.l1.second); // keep consistent interval borders

    // determine gap between overlaps
    if (!isEmpty(sAlocal.second)) {
        // for 2way interactions, the gap must be covered completely by B
        // for 3way interactions, the gap may only be partly covered by B or not covered at all.
        size_t l2start, l2end;
        if (sBlocal.first > sAlocal.first.second)
            l2start = sBlocal.first;
        else
            l2start = sAlocal.first.second;
        if (sBlocal.second < sAlocal.second.first)
            l2end = sBlocal.second;
        else
            l2end = sAlocal.second.first;
        // NOTE: if gap is not covered the interval start is greater than the end and the range will be considered empty
        if (l2end < l2start)
        	region.l2 = empty;
        else
        	region.l2 = make_pair(l2start, l2end);
    } else
        region.l2 = empty;
    empty = make_pair(region.l2.second, region.l2.second); // keep consistent interval borders

    // determine second overlap of A and B
    region.l3 = intersec(sAlocal.second, sBlocal);
    if (SNPDB::isEmpty(region.l3))
		region.l3 = empty; // keep consistent borders!
	empty = make_pair(region.l3.second, region.l3.second); // keep consistent interval borders

    // determine last interval (part right of last overlap of A and B):
    // in the case of 3way, if the exclude range is only for A and not completely covered by B,
    // we need special treatment
    snprange sAlocaltmp = sAlocal.second;
    if (isEmpty(sAlocal.second) || (sBlocal.second < sAlocal.first.second)) { // no gap in A or gap is not covered by B (not even partially)
        region.l4b = sAlocal.second; // could be empty anyway
        sAlocaltmp = sAlocal.first;
    }
    if (sAlocaltmp.second > sBlocal.second) { // end of A is behind end of B
        if (sBlocal.second < sAlocaltmp.first) // gap between A and B
            region.l4 = sAlocaltmp;
        else
            region.l4 = make_pair(sBlocal.second, sAlocaltmp.second);
        region.l4isA = true;
    } else if (sBlocal.second > sAlocaltmp.second) { // end of A is before end of B
        if (sBlocal.first > sAlocaltmp.second) // gap between A and B
            region.l4 = sBlocal;
        else
            region.l4 = make_pair(sAlocaltmp.second, sBlocal.second);
        region.l4isA = false;
    } else
        region.l4 = empty;
    if (SNPDB::isEmpty(region.l4b))
		region.l4b = make_pair(region.l4.second, region.l4.second); // keep interval borders consistent

    // we don't work on the third set for 3rd order interactions
    region.lc = sClocal;

    if (debug) {
		cout << "Local regions:" << endl;
		cout << "l0:\t" << region.l0.first << "\t- " << region.l0.second << "\tis " << (region.l0isA ? "A" : "B") << endl;
		if (order > 2)
			cout << "l0b:\t" << region.l0b.first << "\t- " << region.l0b.second << endl;
		cout << "l1:\t" << region.l1.first << "\t- " << region.l1.second << endl;
		cout << "l2:\t" << region.l2.first << "\t- " << region.l2.second << endl;
		cout << "l3:\t" << region.l3.first << "\t- " << region.l3.second << endl;
		cout << "l4:\t" << region.l4.first << "\t- " << region.l4.second << "\tis " << (region.l4isA ? "A" : "B") << endl;
		if (order > 2) {
			cout << "l4b:\t" << region.l4b.first << "\t- " << region.l4b.second << endl;
			cout << "lc:\t" << region.lc.first << "\t- " << region.lc.second << endl;
		}
//		cout << "--> case " << (region.l0isA ? (region.l4isA ? "a" : "b") : (region.l4isA ? "c" : "d")) << endl;
    }

}

void SNPDB::finishSNP() {
    current_init_snp++;

    current_case_offset.bufptr = &buffer[current_init_snp * snp_size];
    current_case_offset.bit  = 0;

    current_ctrl_offset.bufptr = current_case_offset.bufptr + case_size;
    current_ctrl_offset.bit = 0;
}

snprange SNPDB::parseSetRange(const string &rset_arg) {
    snprange range;
    string rset(rset_arg);

    // maximum range for string "all"
    if (!rset.compare("all")) {
        return make_pair(0, num_snps_global);
    }

    // empty range for string "none"
    if (!rset.compare("none")) {
        return make_pair(0,0);
    }

    // check if string is chromosome-based, i.e. "chrN:a-b" or "N:a-b"
    string chrstr;
    size_t chrpos = rset.find(':');
    if (chrpos != string::npos) {
        chrstr.assign(rset.substr(0,chrpos));
        rset.assign(rset.substr(chrpos+1));
    }

    size_t pos = rset.find('-');
    if (pos == string::npos) { // no range, single element
        size_t idx;
        try {
            idx = stoull(rset);
        } catch (exception &e) {
            StatusFile::addError("No valid range argument.");
            exit(EXIT_FAILURE);
        }
        if (!chrstr.empty()) { // provided chrN:x
            // find SNP in DB
            size_t bp = idx;
            if (!findSNP(chrstr, bp, 0, idx)) {
                stringstream s;
                s << "No SNP found at provided position " << chrstr << ":" << bp;
                StatusFile::addWarning(s.str());
                // setting an empty range
                range.first = idx;
                range.second = idx;
            } else { // range includes only the selected SNP
                range.first = idx;
                range.second = idx+1;
            }
        } else {
            if (idx > num_snps_global) {
                stringstream s;
                s << "Provided index " << idx << " exceeds number of SNPs.";
                StatusFile::addWarning(s.str());
                // setting an empty range
                range.first = num_snps_global;
                range.second = num_snps_global;
            } else {
                range.first = idx-1; // indices are provided 1-based but stored zero-based
                range.second = idx;  // second is excluding
            }
        }
    } else { // range
        size_t idx1;
        size_t idx2;
        bool openrange = false;
        try {
            idx1 = stoull(rset.substr(0,pos)); // parsing will stop at the dash
            if (pos == rset.size()-1) {
                idx2 = num_snps_global; // open range "x-".
                openrange = true;
            }
        } catch (exception &e) {
            StatusFile::addError("No valid range argument.");
            exit(EXIT_FAILURE);
        }
        if (!openrange && chrstr.empty()) {
        	// parse second index
        	try {
        		idx2 = stoull(rset.substr(pos+1)); // parsing after the dash
        	} catch (exception &e) {
        		StatusFile::addError("No valid range argument.");
				exit(EXIT_FAILURE);
        	}
			if (idx2 < idx1) {
				StatusFile::addError(string("No valid range: ") + to_string(idx1) + " - " + to_string(idx2) + ".");
				exit(EXIT_FAILURE);
			}
        }
        if (!chrstr.empty()) { // provided chrN:x-y or chrN:x-chrM:y
            // find SNPs in DB
            size_t bp = idx1;
            findSNP(chrstr, bp, 0, idx1);
            if (idx1 >= num_snps_global) { // chr not found or bp exceeds all SNPs
                stringstream s;
                s << "Chromosome " << chrstr << " not found or position " << bp << " exceeds available SNPs. Set is empty.";
                StatusFile::addWarning(s.str());
                // setting an empty range
                range.first = num_snps_global;
                range.second = num_snps_global;
            } else {
                range.first = idx1; // is including
                if (openrange)
                    range.second = num_snps_global;
                else { // parse next chromosome and index
                	chrpos = rset.find(':',pos+1); // parse after dash
                	if (chrpos != string::npos) { // found a new chromosome string
                		chrstr.assign(rset.substr(pos+1,chrpos-pos-1));
						pos = chrpos;
                	}
                	// parse second bp
                	try {
						idx2 = stoull(rset.substr(pos+1)); // parsing after the dash or double colon
					} catch (exception &e) {
						StatusFile::addError("No valid range argument.");
						exit(EXIT_FAILURE);
					}
                    bp = idx2;
                    if (findSNP(chrstr, bp, idx1, idx2)) // found exactly
                        range.second = idx2+1; // has to be exluding
                    else
                        range.second = idx2; // is excluding anyway
                }
            }
        } else { // provided index range
            if (idx1 > num_snps_global) {
                stringstream s;
                s << "Provided index " << idx1 << " exceeds number of SNPs. Set is empty.";
                StatusFile::addWarning(s.str());
                // setting an empty range
                range.first = num_snps_global;
                range.second = num_snps_global;
            } else {
                if (idx2 > num_snps_global) {
                    stringstream s;
                    s << "Second range argument index " << idx2 << " exceeds number of SNPs.";
                    StatusFile::addWarning(s.str());
                }
                range.first = idx1-1; // indices are provided 1-based but stored zero-based
                range.second = idx2 > num_snps_global ? num_snps_global : idx2;  // second is excluding
            }
        }

    }

    return range;
}

bool SNPDB::findSNP(const string &chr, size_t bp, size_t startidx, size_t &idx) {
    bool foundchr = false;
    for (size_t i = startidx; i < num_snps_global; i++) {
        bool gotchr = snp_info[i].chromosome.compare(chr) == 0;
        if (gotchr)
            foundchr = true;
        if ((gotchr && snp_info[i].pos_bp >= bp) || (!gotchr && foundchr)) { // exceeded bp on chromosome or left chromosome
            idx = i;
            return snp_info[i].pos_bp == bp;
        }
    }
    idx = num_snps_global;
    return false;
}

/*static*/ bool SNPDB::isEmpty(const snprange &r) {
    return r.second <= r.first;
}

/*static*/ size_t SNPDB::lengthOf(const snprange &r) {
    return isEmpty(r) ? 0 : r.second - r.first;
}

/*static*/ bool SNPDB::isLeftOf(const snprange &r1, const snprange &r2) {
    if (isEmpty(r1))
        return false;
    if (isEmpty(r2))
        return true;
    // both ranges are not empty
    return r1.first < r2.first;
}

/*static*/ bool SNPDB::isIncluded(const snprange &r, size_t i) {
    return i >= r.first && i < r.second;
}

/*static*/ snprange SNPDB::intersec(const snprange &r1, const snprange &r2) {
    // NOTE: this works also without checking if an overlap exists first
    snprange r;
    r.first = max(r1.first, r2.first);
    r.second = min(r1.second, r2.second);
    // probably we can skip this part even if we have an empty range anyway...
    if (isEmpty(r)) {
        r.first = 0;
        r.second = 0;
    }
    return r;
}

/*static*/ pair<snprange,snprange> SNPDB::cut(const snprange &r1, const snprange &r2) {
    snprange ret1, ret2;
    snprange tmp = intersec(r1, r2);
    if (isEmpty(tmp)) { // no need to cut
        ret1 = r1;
        ret2.first = 0;
        ret2.second = 0;
    } else { // simply cut the intersection (easy as it is a subset of r1)
        // part left of intersection (could be empty)
        ret1.first = r1.first;
        ret1.second =  tmp.first;
        // part right of intersection (could be empty)
        ret2.first = tmp.second;
        ret2.second = r1.second;
    }
    return make_pair(ret1, ret2);
}

/*static*/ pair<snprange,snprange> SNPDB::unite(const snprange &r1, const snprange &r2) {
    snprange ret1, ret2;

    // sort
    if (!isEmpty(r1) && r1.first <= r2.first) {
        ret1 = r1;
        ret2 = r2;
    } else {
        ret1 = r2;
        ret2 = r1;
    }

    // merge if overlapping
    if (ret2.first <= ret1.second) { // overlap! -> merge
        if (ret2.second > ret1.second)
            ret1.second = ret2.second;
        ret2.first = 0;
        ret2.second = 0;
    } // else: no overlap -> return

    return make_pair(ret1, ret2);
}

/*static*/ vector<snprange> SNPDB::unite(const vector<snprange> &ranges) {
    vector<snprange> ret = ranges; // copy
    bool changed = ranges.size() > 0; // do not use 1 here as the corner case (exactly one empty range) would not be covered then
    while (changed) {
        changed = false;
        for (size_t i2 = 1; i2 < ret.size(); i2++) {
            size_t i1 = i2-1; // always comparing and uniting neighbors (works as the unification also sorts the ranges and puts empties at the back)
            const snprange &r1 = ret[i1];
            const snprange &r2 = ret[i2];
            pair<snprange,snprange> tmp = SNPDB::unite(r1,r2);
            if (tmp.first != r1 || tmp.second != r2) {
                ret[i1] = tmp.first;
                ret[i2] = tmp.second;
                changed = true;
            }
        }
        while (!ret.empty() && SNPDB::isEmpty(ret.back()))
            ret.pop_back();
    }
    return ret;
}

/*static*/ void SNPDB::sort(vector<snprange> &ranges) {
    std::sort(ranges.begin(), ranges.end(), isLeftOf);
    // empty ranges are at the end -> remove
    size_t s = 0;
    for (const auto &r : ranges) {
        if (isEmpty(r))
            break;
        s++;
    }
    ranges.resize(s);
}

/*static*/ vector<snprange> SNPDB::extract(const vector<snprange> &sets, size_t start, size_t length) {
	vector<snprange> extraction;
	if (sets.empty())
		return extraction;
	// find start
	size_t s = 0;
	size_t curr_pos = 0;
	size_t curr_len = 0;
	do {
		curr_pos += lengthOf(sets[s++]);
	} while (s < sets.size() && curr_pos <= start);
	if (curr_pos > start) { // found the start position in the last set
		size_t l = sets[s-1].second - (curr_pos-start) ;
		size_t r = min(l+length, sets[s-1].second);
		extraction.push_back(make_pair(l,r));
		curr_len += r-l;
	}
	// continue until length is reached
	while(s < sets.size() && curr_len < length) {
		extraction.push_back(sets[s]);
		curr_len += lengthOf(sets[s++]);
		if (curr_len > length) { // oh, too much
			extraction.back().second -= curr_len - length;
			curr_len = length; // we're finished
		}
	}
	return extraction;
}

/*static*/ size_t SNPDB::getPairsCnt(size_t a0, size_t a1, size_t b0, size_t b1) {
    // ensure a0 <= b0
    if (b0 < a0) {
        swap(a0, b0);
        swap(a1, b1);
    }
    // the number of pairs is the total rectangle area minus the square part of the overlap plus the "triangle" induced by the overlap
    size_t pairscnt = (a1 - a0) * (b1 - b0); // total area
    pairscnt -= a1 <= b0 ? 0 : ((min(a1, b1) - b0 + 1) * (min(a1, b1) - b0))/2;
    return pairscnt;
}

/*static*/ size_t SNPDB::getPairsCnt(const snprange &r1, const snprange &r2) {
    size_t a0 = r1.first;
    size_t a1 = r1.second;
    size_t b0 = r2.first;
    size_t b1 = r2.second;
    return getPairsCnt(a0, a1, b0, b1);
};

/*static*/ size_t SNPDB::getTriplesCnt(const snprange &r1, const snprange &r2, const snprange &r3) {
    size_t a0 = r1.first;
    size_t a1 = r1.second;
    size_t b0 = r2.first;
    size_t b1 = r2.second;
    size_t c0 = r3.first;
    size_t c1 = r3.second;
    // ensure a0 <= b0 <= c0
    if (b0 < a0) {
        swap(a0, b0);
        swap(a1, b1);
    }
    if (c0 < b0) {
        swap(b0, c0);
        swap(b1, c1);
    }
    if (b0 < a0) {
        swap(a0, b0);
        swap(a1, b1);
    }

    size_t triplescnt = (min(a1, b0) - a0) * getPairsCnt(b0, b1, c0, c1); // indepA x rest
    if (a1 <= b0)
        return triplescnt; // A has no overlaps
    if (b1 < a1) // w.l.o.g. swap rest of A and B such that a1 <= b1
        swap(a1, b1);
    size_t m = min(a1, c0);
    triplescnt += (m - b0) * getPairsCnt(m, b1, c0, c1); // (overlapAB wo. C) x ((rest of B) x C)
    triplescnt += getPairsCnt(b0, m, b0, m) * (c1 - c0); // (pairwise combs. in overlapAB wo. C) x C
    if (a1 <= c0)
        return triplescnt; // no triple overlap
    // w.l.o.g. sort (a1,b1,c1) with already a1<=b1
    if (c1 < b1) {
        swap(b1, c1);
        if (b1 < a1)
            swap(a1, b1);
    }
    triplescnt += (a1 - c0) * getPairsCnt(a1, b1, a1, c1); // overlapABC x ((rest of B) x (rest of C))
    triplescnt += getPairsCnt(c0, a1, c0, a1) * (c1 - a1); // (pairwise combs. in overlapABC) x (rest of C)
    triplescnt += ((a1 - c0) * (a1 - c0 - 1) * (a1 - c0 - 2)) / 6; // triple combs. in overlapABC
    return triplescnt;
}

/*static*/ string SNPDB::printSNPRange(const snprange &r) {
	stringstream s;
	s << r.first << " - " << r.second;
	return s.str();
}
