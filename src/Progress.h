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

#ifndef PROGRESS_H_
#define PROGRESS_H_

#include <time.h>

#include "StatusFile.h"

using namespace std;

class Progress {

public:
    Progress(bool debug_ = false) : debug(debug_) { proc_begin = time(NULL); }

    void start(size_t numtests);
    void updateProgress(size_t numtests_processed);
    void finish(bool ctrlc);

    size_t getNumTests() { return numtests; }
    size_t getNumTestsProcessed() { return numtests_processed; }

private:
    inline void printProgressTime(ostream& out, double progress, time_t begin_timestamp);

    size_t numtests = 0;
    size_t numtests_processed = 0;
    time_t proc_begin; // set time stamp
    double oldprogress = 0.0;

    bool debug;
};

#endif /* PROGRESS_H_ */
