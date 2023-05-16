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

#ifndef GPUENGINE_H
#define GPUENGINE_H

#include <cstdint>
#include <memory>
#include <mutex>
#include <fstream>
#include <sstream>

#include "ResultView.h"
#include "Method.h"

#define GPULOCALDBSIZE 1048576

using namespace std;

class SNPDB;

class GPUEngine

{
public:

	// Eclipse helper
#ifndef RESULT_COMPONENTS
#define RESULT_COMPONENTS 0
#endif

    using ResultViewType = ResultView<>;
    using score_type = ResultViewType::score_type;

    GPUEngine() = delete;
    GPUEngine(const SNPDB &snpdb, const vector<Method> &methods, int index, size_t tableSize, size_t tableBufferSize, ResultView<> view, bool useFPGA, bool debug);
    ~GPUEngine(){};

    void initialize();
    void free();

    void uploadSNPData();
    void runKernel(char* source, char* results, size_t tables_expected);
    int getCUDAIndex() const { return index; }

    void setDumpStream(ostream &dump) {
        dumpptr = &dump;
    }

private:

    const SNPDB &snpdb;
    const vector<Method> &methods;
    unsigned order;
    int index;
    size_t tableSize;
    size_t tableBufferSize;
    size_t resultBufferSize;
    unsigned long dbOffset;

    char *devTables;
    char *devResults;
    ResultView<> resultView;

    bool useFPGA;

    bool debug;

    ostream *dumpptr = 0;
    static mutex dumpmutex;
};

#endif // GPUENGINE_H
