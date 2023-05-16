/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
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

#ifndef PLINKPARSER_H
#define PLINKPARSER_H

#include <vector>
#include <string>

#include "SNPDB.h"

using namespace std;

class PlinkParser
{
public:

    explicit PlinkParser(const string &bim, const string &bed, const string &fam);

    void parse(SNPDB &target, size_t snp_word_size = 1);

    static const size_t input_buffer_size = 1024*1024;
private:

    void parseFam(vector<SNPDB::SampleInfo> &samples, vector<bool> &ignored, size_t &numcases);
    void parseBim(vector<SNPDB::SNPInfo> &snps);
    void parseBed(SNPDB &target, const vector<bool> &ignored_samples);

    string bimfile;
    string bedfile;
    string famfile;

};

#endif // PLINKPARSER_H
