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

#include <iostream>
#include <iomanip>
#include <utility>
#include <chrono>
#include <math.h>
#include <unistd.h>

#include "Progress.h"

void Progress::start(size_t numtests_) {
    numtests = numtests_;
    proc_begin = time(NULL);
    oldprogress = -1.0; // to force a printout of 0%
    updateProgress(0);
}

void Progress::updateProgress(size_t numtests_processed_) {
    numtests_processed = numtests_processed_;
    double progress = numtests_processed_ / (double) numtests;
    // print progress
    double diff = progress - oldprogress;
    if (diff >= 0.01 && (debug || !isatty(STDOUT_FILENO))) {
        cout << numtests_processed << "/" << numtests << " ";
        printProgressTime(cout, progress, proc_begin);
        cout << endl;
        StatusFile::updateStatus(progress);
        oldprogress = progress;
    } else if (diff >= 0.001 && !debug && isatty(STDOUT_FILENO)) {
        cout << "\r                                                                                                  \r";
        printProgressTime(cout, progress, proc_begin);
        StatusFile::updateStatus(progress);
        oldprogress = progress;
    }
}

void Progress::finish(bool ctrlc) {
    oldprogress = -1.0; // to force a printout
    updateProgress(ctrlc ? numtests_processed : numtests);
    if (ctrlc)
        StatusFile::updateStatus(oldprogress, "Terminated");
    else
        StatusFile::updateStatus(1,"Computation finished");
}

void Progress::printProgressTime(ostream& out, double progress, time_t begin_timestamp) {

    out << setprecision(3) << fixed << (progress*100) << " % " << defaultfloat;

    if (progress <= 0.0) {
        out << flush;
        return;
    }

    double totsecs = (difftime(time(NULL), begin_timestamp) / progress);
    double remsecs = totsecs * (1 - progress);

    chrono::seconds ts = chrono::seconds(static_cast<long long int>(round(totsecs)));
    chrono::minutes tm(chrono::duration_cast<chrono::minutes>(ts)); // truncated
    ts -= tm;
    chrono::hours   th(chrono::duration_cast<chrono::hours>(tm)); // truncated
    tm -= th;
    chrono::duration<int,ratio<60*60*24>> td(chrono::duration_cast<chrono::duration<int,ratio<60*60*24>>>(th)); // days, truncated
    th -= td;

    chrono::seconds rs = chrono::seconds(static_cast<long long int>(round(remsecs)));
    chrono::minutes rm(chrono::duration_cast<chrono::minutes>(rs)); // truncated
    rs -= rm;
    chrono::hours   rh(chrono::duration_cast<chrono::hours>(rm)); // truncated
    rm -= rh;
    chrono::duration<int,ratio<60*60*24>> rd(chrono::duration_cast<chrono::duration<int,ratio<60*60*24>>>(rh)); // days, truncated
    rh -= rd;

    out << "\t Rem.time: ";
    if (rd.count())
        out << rd.count() << " days ";
    if (rd.count() || rh.count())
        out << rh.count() << " h ";
    if (rd.count() || rh.count() || rm.count())
        out << rm.count() << " m ";
    out << rs.count() << " s ";

    out << "\t Tot.time: ";
    if (td.count())
        out << td.count() << " days ";
    if (td.count() || th.count())
        out << th.count() << " h ";
    if (td.count() || th.count() || tm.count())
        out << tm.count() << " m ";
    out << ts.count() << " s " << flush;
}
