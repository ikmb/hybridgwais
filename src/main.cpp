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
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <unordered_set>
#include <omp.h>
#include <csignal>
#include <tbb/concurrent_queue.h>

#ifdef USE_CUDA_GPU
#include <cuda_runtime.h>
#endif

#include "hybridsys/Hybridsys.h"
#include "hybridsys/ThreadUtils.h"

#include "version.h"
#include "Args.h"
#include "StatusFile.h"
#include "FPGAConfigurationGWAIS.h"
#include "ThreadPool.h"
#include "FPGAHandler.h"
#include "GPUHandler.h"
#include "GPUIDCreator.h"
#include "CPUProcessor.h"
#include "GPUKernels.h"
#include "ResultHandler.h"
#include "ResultView.h"
#include "SNPDB.h"
#include "PlinkParser.h"
#include "utils.h"
#include "Method.h"
#include "Progress.h"

using namespace std;
using namespace placeholders;

// this is required for a clean shutdown when using FPGAs
static bool term_request = false;
static bool term_wait = false;

StatusFile StatusFile::instance;

/**
 * @brief Signal handler registered for SIGINT reception
 *
 * This function is registered to be called by the operating system if a SIGINT
 * is received (i.e. Ctrl+C). We need to catch this and do a clean process
 * termination when using FPGAs.
 *
 * @param signal The actual signal that has been received
 */
static void signalHandler(int signal) {
    if(signal == SIGINT || signal == SIGTERM) {
        cerr << "\nReceived interrupt. ";
        if(!term_request && term_wait) {
            term_request = true;
            cerr << "Trying to shutdown gracefully, please be patient. Hit Ctrl+C again to terminate immediately." << endl;
        } else {
            cerr << "Terminating..." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

static void verifyConfiguration(hybridsys::Hybridsys& hysys, unsigned order, bool debug) {

    bool useFPGA = hysys.getFPGAs().size() > 0;
    bool useGPU = hysys.getGPUs().size() > 0;

    if (useFPGA && !useGPU) { // FPGA enabled but no GPU
        StatusFile::addError("FPGA acceleration is currently only available together with GPUs.");
        exit(EXIT_FAILURE);
    }

    // get all FPGA configuration blocks
    vector<FPGAConfigurationGWAIS> configs;
    configs.resize(hysys.getFPGAs().size());

    if (useFPGA) {

        for(unsigned i = 0; i < configs.size(); i++) {
            hybridsys::FPGA &f = hysys.getFPGA(i);
            hysys.getFPGAs()[i].createThreadHandle();
            configs[i].parse(f.getConfigurationData());
            hysys.getFPGAs()[i].destroyThreadHandle();
        }

        // check version
        for (auto &cfg: configs) {
            if (get<0>(cfg.getVersion()) != FPGAConfigurationGWAIS::major_required) {
                stringstream ss;
                ss << "Found FPGA version " << get<0>(cfg.getVersion()) << "." << get<1>(cfg.getVersion()) << "rev" << get<2>(cfg.getVersion())
                   << " which is incompatible with this application. Requiring major version " << FPGAConfigurationGWAIS::major_required << ".";
                StatusFile::addError(ss.str());
                exit(EXIT_FAILURE);
            }
            if (get<1>(cfg.getVersion()) < FPGAConfigurationGWAIS::minor_required) {
                stringstream ss;
                ss << "Found FPGA version " << get<0>(cfg.getVersion()) << "." << get<1>(cfg.getVersion()) << "rev" << get<2>(cfg.getVersion())
                   << " which is incompatible with this application. Requiring at least version "
                   << FPGAConfigurationGWAIS::major_required << "." << FPGAConfigurationGWAIS::minor_required << ".";
                StatusFile::addError(ss.str());
                exit(EXIT_FAILURE);
            }
        }

        // assert equality of all FPGA configurations
        if(!all_of(begin(configs) + 1, end(configs), [ &configs ](const FPGAConfigurationGWAIS& other) -> bool { return configs[0] == other; })) {
            StatusFile::addError("All used FPGAs must have the same configurations.");
            exit(EXIT_FAILURE);
        }

        // verify that the requested kernel methods support the underlying FPGA configuration
        switch(order) {
        case 2:
            if(configs[0].getAppID() != FPGAConfigurationGWAIS::APP_2WAY) {
                StatusFile::addError(string("The requested methods are incompatible with the FPGA configuration. Got application ID: ") + to_string(configs[0].getAppID()));
                exit(EXIT_FAILURE);
            }
            break;
        case 3:
            if(configs[0].getAppID() != FPGAConfigurationGWAIS::APP_3WAY) {
                StatusFile::addError(string("The requested methods are incompatible with the FPGA configuration. Got application ID: ") + to_string(configs[0].getAppID()));
                exit(EXIT_FAILURE);
            }
            break;
        default:
            StatusFile::addError("Unsupported method.");
            exit(EXIT_FAILURE);
        }

    } // END if (useFPGA)

    if(debug) {
        if (!useFPGA && !useGPU) {
            cout << "Using CPU-only computation." << endl;
            return;
        }

        if (!useFPGA)
            cout << "Using no FPGAs";
        else {
            cout << "Using " << hysys.getFPGAs().size() << " FPGAs (indices:";
            for(const hybridsys::FPGA& f: hysys.getFPGAs())
                cout << " " << f.getIndex();
            cout << ")";
        }
        cout << " and " << hysys.getGPUs().size() << " GPUs (indices:";
        for(const hybridsys::GPU& g: hysys.getGPUs())
            cout << " " << g.getIndex();
        cout << ")." << endl;
#ifdef USE_CUDA_GPU
        cout << "GPU Capabilities:" << endl;
        for(const hybridsys::GPU& g: hysys.getGPUs()) {
            cout << "GPU " << g.getIndex() << ":" << endl;
            // get CUDA device properties
            cudaDeviceProp deviceProp;
            cudaGetDeviceProperties(&deviceProp, g.getIndex());
            cout << "  Name:                  " << deviceProp.name << endl;
            cout << "  Compute Capability:    " << deviceProp.major << "." << deviceProp.minor << endl;
            cout << "  Clock Rate             " << deviceProp.clockRate << endl;
            cout << "  Total Global Mem:      " << deviceProp.totalGlobalMem << endl;
            cout << "  Total Const. Mem:      " << deviceProp.totalConstMem << endl;
            cout << "  Shared Mem per block:  " << deviceProp.sharedMemPerBlock << endl;
            cout << "  Regs. per block:       " << deviceProp.regsPerBlock << endl;
            cout << "  Glob. Mem. clock rate: " << deviceProp.memoryClockRate << endl;
            cout << "  Glob. Mem. bus width:  " << deviceProp.memoryBusWidth << endl;
            cout << "  L2 cache size:         " << deviceProp.l2CacheSize << endl;
            cout << "  Warp Size:             " << deviceProp.warpSize << endl;
            cout << "  Max threads per block: " << deviceProp.maxThreadsPerBlock << endl;
            cout << "  Max block dimensions:  " << deviceProp.maxThreadsDim[0] << " " << deviceProp.maxThreadsDim[1] << " " << deviceProp.maxThreadsDim[2] << endl;
            cout << "  Max grid dimensions:   " << deviceProp.maxGridSize[0] << " " << deviceProp.maxGridSize[1] << " " << deviceProp.maxGridSize[2] << endl;
            cout << "  Multiprocessor count:  " << deviceProp.multiProcessorCount << endl;
            cout << "  Max threads per MP:    " << deviceProp.maxThreadsPerMultiProcessor << endl;
            cout << "  Shared Mem per MP:     " << deviceProp.sharedMemPerMultiprocessor << endl;
            cout << "  Regs. per MP:          " << deviceProp.regsPerMultiprocessor << endl;
            cout << "  Max 1D texture size:   " << deviceProp.maxTexture1D << endl;
            cout << "  Max 1D surface size:   " << deviceProp.maxSurface1D << endl;
            cout << "  Can map host memory:   " << (deviceProp.canMapHostMemory ? "yes" : "no") << endl;

        }
#endif
        for(unsigned i = 0; i < configs.size(); i++) {
        	cout << "FPGA " << hysys.getFPGA(i).getIndex() << " Configuration: " << endl;
        	FPGAConfigurationGWAIS c = configs[i];
        	cout << c;
        }
    }

}

static void verifyData(const FPGAConfigurationGWAIS &conf, const SNPDB &db, bool useFPGA, bool useGPU, bool debug) {

    stringstream errmsg;
    unsigned long cases = db.getCaseCount();
    unsigned long controls = db.getControlCount();

    if (cases < conf.getMinimumSamples()) {
        stringstream tmp;
        tmp << "<p class='pinfo'>INFO: For best efficiency FPGA design requires at least " << conf.getMinimumSamples() << " cases (got " << cases << ").</p>";
        StatusFile::addInfo(tmp.str());
    }
    if (controls < conf.getMinimumSamples()) {
        stringstream tmp;
        tmp << "<p class='pinfo'>INFO: For best efficiency FPGA design requires at least " << conf.getMinimumSamples() << " controls (got " << controls << ").</p>";
        StatusFile::addInfo(tmp.str());
    }

    cases = db.getCaseCountPadded();
    controls = db.getControlCountPadded();
    unsigned long samples = db.getSampleCountPadded();

    if (cases > conf.getMaximumSamples() && conf.getMaximumSamples() != 0) {
        errmsg << "FPGA design supports at most " << conf.getMaximumSamples() << " cases or controls (got " << cases << ", including padding).";
        StatusFile::addError(errmsg.str());
        exit(EXIT_FAILURE);
    }
    if (controls > conf.getMaximumSamples() && conf.getMaximumSamples() != 0) {
        errmsg << "FPGA design supports at most " << conf.getMaximumSamples() << " cases or controls (got " << controls << ", including padding).";
        StatusFile::addError(errmsg.str());
        exit(EXIT_FAILURE);
    }

    if(db.getSNPSize() * 2 > GPULOCALDBSIZE && !useFPGA && useGPU) {
        errmsg << "GPU only run supports at most " << (GPULOCALDBSIZE*2) << " samples in total (got " << db.getSampleCountPadded() << ", including padding).";
        StatusFile::addError(errmsg.str());
        exit(EXIT_FAILURE);
    }

    if(debug) {
        cout << "*Padded* data set consists of " << cases << " cases and " << controls << " controls (" << samples << " samples in total)." << endl;
    }
}

static size_t getGPUResultSize(const vector<Method>& methods) {
    size_t gpu_result_size = methods[0].getOrder() * sizeof(ResultView<>::id_type);
    // TODO why is the alignment already required here?
    if(gpu_result_size % ResultView<>::getAlignment() != 0)
        gpu_result_size += (ResultView<>::getAlignment() - (gpu_result_size % ResultView<>::getAlignment()));

    gpu_result_size += Method::getScoreFields(methods) * sizeof(ResultView<>::score_type);
    if(gpu_result_size % ResultView<>::getAlignment() != 0)
        gpu_result_size += (ResultView<>::getAlignment() - (gpu_result_size % ResultView<>::getAlignment()));

    return gpu_result_size;
}


int main(int argc, char *argv[]) {
    ThreadUtils::setThreadName("hygwais");

    cout << "--------------------------------------------------------" << endl;
    cout << "--- HybridGWAIS: Fast Exhaustive Interaction Detection" << endl;
    cout << "---" << endl;
    cout << "--- " << prog_version << endl;
    cout << "--- (compiled on " <<  prog_timestamp << ")" << endl;
    cout << "---" << endl;
    cout << "--- Copyright (C) 2018-2023 by Lars Wienbrandt," << endl;
    cout << "--- Institute of Clinical Molecular Biology, Kiel University" << endl;
    cout << "--- Distributed under the GNU GPLv3 license." << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << endl;

    if (argc > 1) {
        cout << "Provided arguments:" << endl;
        cout << "  ";
        size_t cpos = 0;
        for (int a = 1; a < argc; a++) { // skip first argument (call to binary)
            if (a>1) {
                if (argv[a][0] == '-') {
                    cout << endl << "  ";
                    cpos = 0;
                } else {
                    for (; cpos < 18; cpos++) {
                        cout << ' ';
                    }
                }
            }
            size_t p = string(argv[a]).find_last_of('/');
            if (p == string::npos) {
                cout << argv[a] << " ";
                cpos += string(argv[a]).size()+1;
            } else
                cout << "***" << string(argv[a]).substr(p); // mask all paths, do not need to update cpos here
        }
        cout << "\n" << endl;
    } // else: no arguments provided. will exit in Args.

    Args args = Args::parseArgs(argc, argv);
    bool yaml = args.yaml;

    // register Ctrl+C handler for proper cleanup on interruption
    signal(SIGINT, signalHandler);

    {
        ifstream f("/proc/self/status");
        string line;
        int tracerPid = 0;
        static const char token[] = "TracerPid:";

        while(getline(f, line)) {
            size_t pos = line.find(token);
            if(pos != string::npos) {
                tracerPid = atoi(line.data() + pos + sizeof(token));
            }
        }
        if (tracerPid > 0)
        	StatusFile::addWarning("Detected attached debugger!");
    }

    vector<Method> &methods = args.methods;
    if (methods.size() == 0) {
        StatusFile::addError("No method defined.");
        exit(EXIT_FAILURE);
    }
    unsigned order = methods[0].getOrder();
    for (Method &m : methods) {
        if (m.getOrder() != order) {
            StatusFile::addError("Mixed orders! Only methods of the same order are allowed.");
            exit(EXIT_FAILURE);
        }
        m.setDetails(args.decompose);
    }
    // add LD to methods, if desired
    if (args.ld) {
    	if (order == 2)
    		methods.emplace_back(Method::Type::LD);
    	else
    		methods.emplace_back(Method::Type::LD3);
    	methods.back().setDetails(args.decompose);
    }

    // initialize platform
    StatusFile::updateStatus(0,"Initializing");

#ifdef USE_CUDA_GPU
#ifdef USE_AD_FPGA
    hybridsys::Hybridsys hysys(args.fpgas, args.gpus);
#else // GPU only
    hybridsys::Hybridsys hysys(vector<unsigned>(), args.gpus);
#endif
#else // CPU only
    hybridsys::Hybridsys hysys;
#endif

    bool useFPGA = hysys.getFPGAs().size() > 0;
    bool useGPU = hysys.getGPUs().size() > 0;

    if (useFPGA && order > 2) {
    	StatusFile::addError("FPGA usage for 3rd order methods is disabled.");
    	exit(EXIT_FAILURE);
    }

    if (useFPGA) {
        cout << "Using FPGA(s):";
        for (const auto &fpga : hysys.getFPGAs())
            cout << " " << fpga.getIndex();
        cout << endl;
    }

    if (useGPU) {
        cout << "Using GPU(s):";
        for (const auto &gpu : hysys.getGPUs())
            cout << " " << gpu.getIndex();
        cout << endl;
    }

    verifyConfiguration(hysys, order, args.debug);
#if defined USE_CUDA_GPU && defined USE_AD_FPGA
    chrono::milliseconds fpga_timeout(args.timeout);
#else
    chrono::milliseconds fpga_timeout(0);
#endif

    const unsigned num_buffers_per_fpga = useFPGA ? args.nbuffersFPGA : (useGPU ? hysys.getGPUs().size()+1 : 0); // GPU-only: used for ID creation
    const unsigned num_buffers_per_gpu = useGPU ? args.nbuffersGPU : 0;
	if (useFPGA && num_buffers_per_fpga == 0) {
		StatusFile::addError("Need to specify at least one transfer buffer for FPGA->GPU!");
		exit(EXIT_FAILURE);
	}
	if (useGPU && num_buffers_per_gpu == 0) {
		StatusFile::addError("Need to specify at least one transfer buffer for GPU->ResultProcessing!");
		exit(EXIT_FAILURE);
	}

	size_t excluderange = args.excluderange;

	if (order == 2) {
	    if (args.setC_str.compare("all"))
	        StatusFile::addWarning("Ignoring provided third user region as no third order method was selected.");
	    args.setC_str = string("none"); // required for 2way methods
	}

	FPGAConfigurationGWAIS conf(useFPGA ? hysys.getFPGA(0).getConfigurationData() : NULL);
    SNPDB db(conf.getMinimumSamples(), conf.getMinimumSamples(), args.debug);
    db.setUserRegions(args.setA_str, args.setB_str, args.setC_str, args.exset_str, args.useexsetA, args.excluderange, order);

    string bim = args.bim;
    string bed = args.bed;
    string fam = args.fam;
    size_t p;
    cout << "\nFiles:" << endl;
    p = bed.find_last_of('/');
    cout << " bed: " << (p == string::npos ? bed : string("***")+bed.substr(p)) << endl;
    p = bim.find_last_of('/');
    cout << " bim: " << (p == string::npos ? bim : string("***")+bim.substr(p)) << endl;
    p = fam.find_last_of('/');
    cout << " fam: " << (p == string::npos ? fam : string("***")+fam.substr(p)) << endl;

    // load data
    cout << "\nLoading data... " << flush;
	StatusFile::updateStatus(0,"Loading data");

    PlinkParser parser(bim, bed, fam);
    if (useFPGA)
        parser.parse(db, conf.getSnpWordSize());
    else if (useGPU)
        parser.parse(db, 512); // TODO the word size should be set according to GPU memory specifications (memory bus width in bytes or so)
    else // host-only
        parser.parse(db);

    cout << "done.\n" << endl;

    // prepare status
    StatusFile::addInfo("<h3>General Information:</h3>", !yaml);
    stringstream infostr;
    infostr << "<table><tr><td><b> SNPs (all):  </b></td><td><b>" << db.getGlobalSNPCount() << "</b></td></tr>\n";
    infostr <<        "<tr><td><b> SNPs (used): </b></td><td><b>" << db.getLocalSNPCount() << "</b></td></tr>\n";
    infostr <<        "<tr><td><b> Samples:   </b></td><td><b>" << db.getSampleCount() << "</b></td></tr>\n";
    infostr <<           "<tr><td>  Cases:    </td><td>" << db.getCaseCount() << "</td></tr>\n";
    infostr <<           "<tr><td>  Controls: </td><td>" << db.getControlCount() << "</td></tr></table>\n";
    StatusFile::addInfo(infostr.str(), !yaml);

    if (yaml) {
        stringstream istr;
        istr << "\n"
                << "  SNPs:\n"
				<< "    global: " << db.getGlobalSNPCount() << "\n"
                << "    local:  " << db.getLocalSNPCount() << "\n"
                << "  samples:\n"
                << "    total:    " << db.getSampleCount() << "\n"
                << "    cases:    " << db.getCaseCount() << "\n"
                << "    controls: " << db.getControlCount() << "\n";
        StatusFile::addInfoYAML("General information", istr.str());
    }

    verifyData(conf, db, useFPGA, useGPU, args.debug);
    if (args.debug)
        cout << "SNP size: " << db.getSNPSize() << endl;

    // output filename base -> prefix is built from first method to process
    stringstream outfiless;
    outfiless << args.outprefix;
    for (const auto &m : methods)
        outfiless << "." << m.getShortName();

    GPUEngine::score_type threshold = args.threshold;
    unsigned long bestresults = args.nresults;

    // User regions
    SNPDB::RegionInfo regioninfo = db.getRegionInfo();
    snprange usersetAg = regioninfo.setA_g;
    snprange usersetBg = regioninfo.setB_g;
    snprange usersetCg = regioninfo.setC_g;
    {
        stringstream setstream;
        setstream << "<h3>SNP sets:</h3>\n";
        setstream << "<table><tr><td> Set A: </td>";
        if (usersetAg.first >= usersetAg.second)
            setstream << "<td>empty</td></tr>\n";
        else
            setstream << "<td>" << (usersetAg.first + 1)
                << "</td> <td>-</td> <td>" << usersetAg.second
                << "</td>\t<td>==</td> <td>" << db.getSNPInfo(usersetAg.first).chromosome << ":" << db.getSNPInfo(usersetAg.first).pos_bp
                << "</td> <td>-</td> <td>" << db.getSNPInfo(usersetAg.second-1).chromosome << ":" << db.getSNPInfo(usersetAg.second-1).pos_bp << "</td></tr>\n"; // 1-based indices and end is included
        setstream << "<tr><td> Set B: </td>";
        if (usersetBg.first >= usersetBg.second)
            setstream << "<td>empty</td></tr>\n";
        else
            setstream << "<td>" << (usersetBg.first + 1)
                << "</td> <td>-</td> <td>" << usersetBg.second
                << "</td>\t<td>==</td> <td>" << db.getSNPInfo(usersetBg.first).chromosome << ":" << db.getSNPInfo(usersetBg.first).pos_bp
                << "</td> <td>-</td> <td>" << db.getSNPInfo(usersetBg.second-1).chromosome << ":" << db.getSNPInfo(usersetBg.second-1).pos_bp << "</td></tr>\n"; // 1-based indices and end is included
        if (order == 3) {
            setstream << "<tr><td> Set C: </td>";
            if (usersetCg.first >= usersetCg.second)
                setstream << "<td>empty</td></tr>\n";
            else
                setstream << "<td>" << (usersetCg.first + 1)
                    << "</td> <td>-</td> <td>" << usersetCg.second
                    << "</td>\t<td>==</td> <td>" << db.getSNPInfo(usersetCg.first).chromosome << ":" << db.getSNPInfo(usersetCg.first).pos_bp
                    << "</td> <td>-</td> <td>" << db.getSNPInfo(usersetCg.second-1).chromosome << ":" << db.getSNPInfo(usersetCg.second-1).pos_bp << "</td></tr>\n"; // 1-based indices and end is included
        }
        setstream << "</table>";
        StatusFile::addInfo(setstream.str(), !yaml);

        if (yaml) {
            stringstream sstr;
            sstr << "\n"
                      << "    setA: ";
            if (usersetAg.first >= usersetAg.second)
                sstr << "empty\n";
            else {
                sstr << "\n"
                      << "      beginIdx: " << (usersetAg.first + 1) << "\n"
                      << "      endIdx: " << usersetAg.second << "\n"
                      << "      beginGenPos: " << db.getSNPInfo(usersetAg.first).chromosome << ":" << db.getSNPInfo(usersetAg.first).pos_bp << "\n"
                      << "      endGenPos: " << db.getSNPInfo(usersetAg.second-1).chromosome << ":" << db.getSNPInfo(usersetAg.second-1).pos_bp << "\n";
            }
            sstr << "    setB: ";
            if (usersetBg.first >= usersetBg.second)
                sstr << "empty\n";
            else {
                sstr << "\n"
                      << "      beginIdx: " << (usersetBg.first + 1) << "\n"
                      << "      endIdx: " << usersetBg.second << "\n"
                      << "      beginGenPos: " << db.getSNPInfo(usersetBg.first).chromosome << ":" << db.getSNPInfo(usersetBg.first).pos_bp << "\n"
                      << "      endGenPos: " << db.getSNPInfo(usersetBg.second-1).chromosome << ":" << db.getSNPInfo(usersetBg.second-1).pos_bp << "\n";
            }
            if (order == 3) {
                sstr << "    setC: ";
                if (usersetCg.first >= usersetCg.second)
                    sstr << "empty\n";
                else {
                    sstr << "\n"
                          << "      beginIdx: " << (usersetCg.first + 1) << "\n"
                          << "      endIdx: " << usersetCg.second << "\n"
                          << "      beginGenPos: " << db.getSNPInfo(usersetCg.first).chromosome << ":" << db.getSNPInfo(usersetCg.first).pos_bp << "\n"
                          << "      endGenPos: " << db.getSNPInfo(usersetCg.second-1).chromosome << ":" << db.getSNPInfo(usersetCg.second-1).pos_bp << "\n";
                }
            }
            StatusFile::addInfoYAML("SNP sets", sstr.str());
        }

        if (SNPDB::isEmpty(usersetAg) || SNPDB::isEmpty(usersetBg) || (order == 3 && SNPDB::isEmpty(usersetCg))) {
            StatusFile::addError("Empty SNP sets.");
            exit(EXIT_FAILURE);
        }
    }

    bool useexset = regioninfo.useexset;
    bool useexsetA = regioninfo.useexsetA;
    snprange userexsetg = regioninfo.exset_g;
    if (useexset || excluderange > 0) { // using user specified SNP ranges or an excluderange in bp
        stringstream sstr;
        sstr << "<h3>User specified SNP exclude sets:</h3>" << endl;
        if (useexset) {
            sstr << "<table><tr><td> Exclude set" << (useexsetA ? " (only for A)" : "") << ": </td>";
            if (userexsetg.first >= userexsetg.second)
                sstr << "<td>empty</td></tr>\n";
            else
                sstr << "<td>" << (userexsetg.first + 1)
                << "</td> <td>-</td> <td>" << (userexsetg.second)
                << "</td>\t<td>==</td> <td>" << db.getSNPInfo(userexsetg.first).chromosome << ":" << db.getSNPInfo(userexsetg.first).pos_bp
                << "</td> <td>-</td> <td>" << db.getSNPInfo(userexsetg.second-1).chromosome << ":" << db.getSNPInfo(userexsetg.second-1).pos_bp << "</td></tr>\n";
            sstr << "</table>";
        }
        if (excluderange > 0)
            sstr << "<table><tr><td> General exclude range: </td><td>" << excluderange << " bp</td></tr></table>\n";

        StatusFile::addInfo(sstr.str(), !yaml);

        if (SNPDB::isEmpty(userexsetg) && useexset)
            StatusFile::addWarning("Exclude set is empty.");
    }

    if (yaml) {
        stringstream sstr;
        sstr << "\n";
		sstr << "    excludeset: ";
        if (userexsetg.first >= userexsetg.second)
            sstr << "empty\n";
        else {
            sstr << "\n"
                 << "      beginIdx: " << (userexsetg.first + 1) << "\n"
                 << "      endIdx: " << userexsetg.second << "\n"
                 << "      beginGenPos: " << db.getSNPInfo(userexsetg.first).chromosome << ":" << db.getSNPInfo(userexsetg.first).pos_bp << "\n"
                 << "      endGenPos: " << db.getSNPInfo(userexsetg.second-1).chromosome << ":" << db.getSNPInfo(userexsetg.second-1).pos_bp << "\n"
                 << "      is only for A: " << (useexsetA ? "true" : "false") << "\n";
        }

        sstr << "    excluderange: " << excluderange << "\n";

        StatusFile::addInfoYAML("SNP exclude sets", sstr.str());
    }

    // buffer sizes and result size
    size_t gpuresultsize = useGPU ? getGPUResultSize(methods) : 0;
    size_t fpgatablesize = useGPU ? (!useFPGA ? (order == 2 ? FPGAHandler::GPUONLY_IDSIZE_2WAY : FPGAHandler::GPUONLY_IDSIZE_3WAY) : conf.getTableSize()) : 0;
    size_t fpgabuffersize = useGPU ? args.buffersize : 0;
    size_t gpubuffersize = useGPU ? (!useFPGA ?
            gpuresultsize * (fpgabuffersize / gpuresultsize) : // take the provided buffer size as GPU buffer size, but floored to exactly fit a multiple of the results
            gpuresultsize * (fpgabuffersize / fpgatablesize))  // the size of the results for each table from the FPGA
            : 0; // not specified for CPU-only run
    if (useGPU && !useFPGA)
        fpgabuffersize = (fpgabuffersize / gpuresultsize) * fpgatablesize; // adjusted to keep only the IDs for each result

    // define a view on the results
    ResultView<> view(order, Method::getScoreFields(methods), gpubuffersize); // gpubuffersize == 0 for CPU-only runs

    ResultHandler<typename ResultView<>::id_type, typename ResultView<>::score_type> resultHandler(
    			db,
                outfiless.str(),
                useGPU ? num_buffers_per_gpu : omp_get_max_threads(), // provide one thread for each available buffer or all available threads for CPU-only runs
                threshold,
                bestresults,
                view,
                methods,
                args
                );

    // main() return value
    int retval = EXIT_FAILURE;
    // progress
    Progress progress(args.debug);

    if (useGPU) {
#ifdef USE_CUDA_GPU
        ///////////////////////////////
        //  hardware accelerated run //
        ///////////////////////////////

        hybridsys::BufferFactory<hybridsys::FPGABuffer> fpgaBufferFactory(fpgabuffersize);
        hybridsys::BufferFactory<hybridsys::CUDABuffer> gpuBufferFactory(gpubuffersize);

        if (args.debug) {
            cout << "FPGA buffer size: " << fpgabuffersize << " CUDA buffer size: " << gpubuffersize
                << " Items per buffer: " << (gpubuffersize / gpuresultsize) << endl;
        }

        if (useFPGA) {
            for(unsigned i = 0; i < hysys.getFPGAs().size(); ++i)
                for(unsigned j = 0; j < num_buffers_per_fpga; ++j) {
                    fpgaBufferFactory.preallocateBuffer();
                }
        } else {
            for(unsigned j = 0; j < num_buffers_per_fpga; ++j) {
                fpgaBufferFactory.preallocateBuffer(); // allocate buffers for the IDs
            }
        }


        for(unsigned i = 0; i < hysys.getGPUs().size(); ++i) {

            cudaSetDevice(hysys.getGPU(i).getIndex());

            for(unsigned j = 0; j < num_buffers_per_gpu; ++j) {
                gpuBufferFactory.preallocateBuffer();
            }
        }

        if (args.debug)
            cout << "Preloaded buffers." << endl;


        // set up transmission queues
        tbb::concurrent_bounded_queue<shared_ptr<hybridsys::FPGABuffer>> dummyQueue, gpuQueue;
        tbb::concurrent_bounded_queue<shared_ptr<hybridsys::CUDABuffer>> resultProcessorQueue;

        // set up queue handlers
        GPUHandler gpuHandler(hysys, db, gpuBufferFactory, fpgatablesize, fpgabuffersize, methods, view, useFPGA, args.debug);
        cout << "Prepared GPU(s)." << endl;

        if (useFPGA) {
            StatusFile::updateStatus(0,"Waiting for FPGA");
            term_wait = true; // need proper shutdown procedure when using FPGAs
        }

        // we need either the FPGAHandler or the GPUIDCreator depending on FPGA usage
        FPGAHandler fpgaHandler(hysys, fpga_timeout, fpgaBufferFactory, db, methods, progress, term_request, args.debug);
        GPUIDCreator gpuidcreator(fpgaBufferFactory, db, methods, progress, args.debug);

        if (useFPGA)
            cout << "Prepared FPGA(s)." << endl;

        // one thread per channel (1 ID creator if we don't want to use the FPGAs -> This is definitely a TODO as this is the bottleneck for GPU-only processing!)
        int fpgaFeederCount = useFPGA ? (hysys.getFPGAs().size() * conf.getNumChains()) : 1;


        // starting the FPGA handler threads (or the ID creator thread if no FPGAs are used)
        ThreadPool<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::FPGABuffer>> fpgaProcessor(fpgaFeederCount,
                                                                                           dummyQueue,
                                                                                           gpuQueue,
                                                                                           useFPGA ? (ThreadPool<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::FPGABuffer>>::threadfunc_type)
                                                                                                       bind(&FPGAHandler::fetchTables, &fpgaHandler, _1, _2, _3, _4)
                                                                                                   : (ThreadPool<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::FPGABuffer>>::threadfunc_type)
                                                                                                       bind(&GPUIDCreator::createIDs, &gpuidcreator, _1, _2, _3, _4));

        // starting GPU handler threads (one thread per GPU)
        ThreadPool<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::CUDABuffer>> gpuProcessor(hysys.getGPUs().size(),
                                                                                           gpuQueue,
                                                                                           resultProcessorQueue,
                                                                                           (ThreadPool<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::CUDABuffer>>::threadfunc_type)
                                                                                               bind(&GPUHandler::process<shared_ptr<hybridsys::FPGABuffer>, shared_ptr<hybridsys::CUDABuffer>>,
                                                                                                    &gpuHandler, _1, _2, _3, _4));

        // starting result handler threads (one thread per GPU buffer)
        ThreadPool<shared_ptr<hybridsys::CUDABuffer>, shared_ptr<hybridsys::FPGABuffer>> resultProcessor(num_buffers_per_gpu,
                                                                                           resultProcessorQueue,
                                                                                           dummyQueue,
                                                                                           (ThreadPool<shared_ptr<hybridsys::CUDABuffer>, shared_ptr<hybridsys::FPGABuffer>>::threadfunc_type)
                                                                                               bind(&ResultHandler<typename ResultViewSpec::id_type, typename ResultViewSpec::score_type>::process<shared_ptr<hybridsys::CUDABuffer>,
                                                                                                       shared_ptr<hybridsys::FPGABuffer>>, &resultHandler, _1, _2, _3, _4));

        cout << "\nRunning..." << endl;
        StatusFile::updateStatus(0,"Running");

        vector<thread> fpgacrtab;
        if (useFPGA)
            fpgacrtab.emplace_back(&FPGAHandler::createTables, std::ref(fpgaHandler));

        if(!isatty(STDOUT_FILENO) || args.debug)
            cout << endl;

        bool terminated = false;
        while(!terminated) {
            this_thread::sleep_for(chrono::milliseconds(500)); // poll twice a second

            if(term_request || fpgaHandler.timeoutOccurred() || fpgaHandler.isFinished() || gpuidcreator.isFinished()) {

                // cancel all DMA transactions on FGPA
                if (term_request || fpgaHandler.timeoutOccurred()) {
                    term_request = true; // a timeout is handled the same way as a user termination
                    for (auto &fpga : hysys.getFPGAs()) {
                        fpga.abort();
//                        fpga.cancel(0); // write channel
//                        for (unsigned int c = 1; c <= conf.getNumChains(); c++) {
//                            fpga.cancel(c); // read channels
//                        }
                    }
                }

                fpgaProcessor.terminate(false);
                if (args.debug && useFPGA)
                    cout << "\nFPGAs have terminated. FPGA handler was aborted: " << (fpgaHandler.isFinished() ? "no" : "yes") << endl;
                gpuProcessor.terminate(fpgaHandler.isFinished() || gpuidcreator.isFinished()); // true on normal finish, false on user interrupt or timeout
                if (args.debug)
                    cout << "\nGPUs have terminated. " << endl;
                resultProcessor.terminate(fpgaHandler.isFinished() || gpuidcreator.isFinished()); // true on normal finish, false on user interrupt or timeout
                if (args.debug)
                    cout << "Result processing has terminated. " << endl;
                if (useFPGA)
                    fpgacrtab[0].join();
                terminated = true;
            }

        }
        term_wait = false; // critical part for user termination is over
#endif
    } else { // !useGPU
#ifndef USE_CUDA_GPU
        (void) blockIdx;
        (void) blockDim;
        (void) threadIdx;
#endif
        ////////////////////////////
        //      CPU-only run      //
        ////////////////////////////

        cout << "\nRunning..." << endl;
        StatusFile::updateStatus(0,"Running");

        CPUProcessor proc(db, methods, progress, view, resultHandler, args.debug);
        proc.process();

    } // end if(useGPU) ... else ...

    progress.finish(term_request);
    if (!term_request) { // finished
        StatusFile::updateStatus(1,"Writing results");
        string rfile = resultHandler.getResultFileName();
        p = rfile.find_last_of('/');
        cout << "\n\nResult file: " << (p == string::npos ? rfile : string("***")+rfile.substr(p)) << endl;
        cout << "Writing results..." << endl;
        string cmdline;
        for (int i = 0; i < argc; i++) { // all args including the executable call
            cmdline += " ";
            p = string(argv[i]).find_last_of('/');
            if (p == string::npos)
                cmdline += argv[i];
            else
                cmdline += string("***") + string(argv[i]).substr(p);
        }
        resultHandler.flush(cmdline);
        StatusFile::updateStatus(1,"Computation finished");
        retval = EXIT_SUCCESS;
    } else { // User terminate or timeout
        cout << "Terminated." << endl;
    }

    return retval;
}
