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
#include <cmath>
#include <thread>

#include <boost/program_options.hpp>

#include "version.h"
#include "Args.h"
#include "GPUEngine.h"
#include "Method.h"
#include "StatusFile.h"

namespace bpo = boost::program_options;

using namespace std;
using namespace bpo;

/**
 * @brief Constructs, parses and verifies all command-line arguments
 *
 * This function constructs the Args object and parses all given command-line
 * arguments. If unknown arguments and/or missing obligatory arguments are
 * detected, this function does not return and instead prints an appropriate
 * help message and calls exit(), executing all exit handlers registered
 * up to the parseArgs call.
 *
 * @param argc argument count
 * @param argv argument vector
 * @return an rvalue reference of a new Args object containing all defined or defaulted CLI arguments
 */
Args Args::parseArgs(int argc, char *argv[]) {
    Args args {argc, argv};

    if(argc <= 1 || args.count("help")) {
        args.printHelp(argv[0], cout);
        exit(EXIT_SUCCESS);
    }

    if(args.count("version")) {
        printVersion();
        exit(EXIT_SUCCESS);
    }

    if(!(args.count("snpfile") || (args.count("bim") && args.count("bed") && args.count("fam")))) {
        StatusFile::addError("You need to specify an input file basename OR individual file names with --{bim,bed,fam}. You may not specify both or neither.");
        exit(EXIT_FAILURE);
    }

    if(!args.count("method")) {
        StatusFile::addError("Missing required option --method/-m");
        exit(EXIT_FAILURE);
    }

    if(!args.count("stat") && args.count("yaml")) {
        StatusFile::addError("You need to specify a status file base with --stat if you want YAML messages.");
        exit(EXIT_FAILURE);
    }

    // set variables
    args.parseVars();

    return args;
}

ostream &operator<<(ostream &out, const Args &args) {
    variables_map::const_iterator it;

    long name_width = 0;
    long value_width = 0;
    long orig_width = out.width();

    // collect field widths
    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        long this_width = static_cast<long>(it->first.length());

        if(this_width > name_width) {
            name_width = this_width;
        }

        this_width = 0;

        if(it->second.value().type() == typeid(string)) {
            this_width = static_cast<long>(it->second.as<string>().length());
        }

        if(it->second.value().type() == typeid(int)) {
            this_width = static_cast<long>(log10(it->second.as<int>()));
        }

        if(this_width > value_width) {
            value_width = this_width;
        }
    }

    // dump option values
    out.setf(ios::adjustfield, ios::left);

    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        out.width(name_width + 2); // more space for the colon
        out << (it->first + ":") << endl;

        out.width(value_width);

        if(it->second.value().type() == typeid(string)) {
            out << it->second.as<string>() << endl;
        } else if(it->second.value().type() == typeid(int)) {
            out << it->second.as<int>() << endl;
        } else {
            out << "(unknown)" << endl;
        }
    }

    resetiosflags(ios::adjustfield);
    out.width(orig_width);

    return out;
}

Args::Args(int argc, char *argv[]) :
    opts_regular("Program options"),
	opts_regular_hybrid("Additional program options for hybrid architecture"),
    opts_hidden("Hidden options (only visible in debug mode)"),
	opts_hidden_hybrid("Additional hidden options for hybrid architecture")
	{

    opts_regular.add_options()
    ("help,h", "produce this help message and terminates")
    ("version", "prints version information and terminates")
    ("output,o", value<string>(&outprefix), "output file(s prefix, if applicable)")
    ("stat", value<string>(&statfile), "file for current status output (instead of stdout)")
    ("method,m", value<vector<Method>>(&methods)->multitoken(), "list of test methods to apply concurrently (see below). It is not allowed to mix 2nd and 3rd order methods!")
	("ld", "Compute linkage disequilibrium (r2) alongside the selected methods. For 3rd-order methods all pairwise LDs are calculated.")
	("ldfilter", value<GPUEngine::score_type>(&ldfilter)->default_value(0.0), "Computes linkage disequilibrium and filters according to the provided values. Only results with an r2 <= <arg> pass the filter. For 3rd-order methods only those triples where all pairwise LDs are less than <arg> pass the filter. This implicitly activates --ld.")
    ("nresults,n", value<size_t>(&nresults)->default_value(100000UL), "select only the N best results")
	("decompose", "decompose result values (e.g. information gain is decomposed into mutual information subvalues)")
    ("threshold,t", value<GPUEngine::score_type>(&threshold)->default_value(numeric_limits<GPUEngine::score_type>::quiet_NaN()), "threshold for result selection, but will be capped at the 1 million best results if not stated otherwise with -n")
    ("bed", value<string>(&bed), "BED file input (can be used instead of the <SNP input file> basename")
    ("bim", value<string>(&bim), "BIM file input (can be used instead of the <SNP input file> basename")
    ("fam", value<string>(&fam), "FAM file input (can be used instead of the <SNP input file> basename")
    ("setA", value<string>(&setA_str)->default_value("all"), "Subset of SNPs for interaction tests at first position in pair or triple. chr:bpStart-[chr:]bpEnd, 1-based SNP index range, or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("setB", value<string>(&setB_str)->default_value("all"), "Subset of SNPs for interaction tests at second position in pair or triple. chr:bpStart-[chr:]bpEnd, 1-based SNP index range, or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("setC", value<string>(&setC_str)->default_value("all"), "Subset of SNPs for interaction tests at third position in triple. chr:bpStart-[chr:]bpEnd, 1-based SNP index range, or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("exset", value<string>(&exset_str)->default_value("none"), "Subset of SNPs for interaction tests to be excluded. chr:bpStart-[chr:]bpEnd, 1-based SNP index range, or \"none\" (default). Range may as well be a single SNP. Cannot be used together with --exsetA")
    ("exsetA", value<string>(&exsetA_str)->default_value("none"), "Subset of SNPs for interaction tests to be excluded only at first position in pair or triple. chr:bpStart-bpEnd, 1-based SNP index range, or \"none\" (default). Range may as well be a single SNP. Cannot be used together with --exset")
//    ("exsetB", value<string>(&exsetB_str)->default_value("none"), "Subset of SNPs for interaction tests to be excluded at second position in pair or triple. chr:bpStart-bpEnd, 1-based SNP index range, or \"none\" (default). Range may as well be a single SNP.")
//    ("exsetC", value<string>(&exsetC_str)->default_value("none"), "Subset of SNPs for interaction tests to be excluded at third position in triple. chr:bpStart-bpEnd, 1-based SNP index range, or \"none\" (default). Range may as well be a single SNP.")
    ("exrange", value<size_t>(&excluderange)->default_value(0), "exclude test if a pair (pairwise tests) or one pair of a triple (3way tests) is within this range (in bp)")
    ;

#ifdef USE_CUDA_GPU
#ifdef USE_AD_FPGA
    opts_regular_hybrid.add_options()
    ("fpga", value<vector<unsigned>>(&fpgas)->multitoken()->default_value(vector<unsigned>(),""), "restrict to these FPGAs (separate indices by whitespace, defaults to no FPGA acceleration)")
    ;
#endif
    opts_regular_hybrid.add_options()
    ("gpu", value<vector<unsigned>>(&gpus)->multitoken()->default_value(vector<unsigned>(),""), "restrict to these GPUs (separate indices by whitespace, defaults to no GPU acceleration)")
    ;
#endif

    opts_hidden.add_options()
    ("debug", "produce lots of debug output")
    ("snpfile", value<string>(&snpfilebase), "SNP input file (filename base for .bim/.bed/.fam (automatic parameter for positional argument #1)")
    ("yaml", "produce info/warning/error files in YAML format")
	("snpindex", "put SNP indices in result file")
	;

#ifdef USE_CUDA_GPU
#ifdef USE_AD_FPGA
    opts_hidden_hybrid.add_options()
    ("buffersize", value<size_t>(&buffersize)->default_value(256*1024*1024ull), "Size for transmission buffers (FPGA->GPU) in bytes.")
    ("buffersFPGA", value<unsigned>(&nbuffersFPGA)->default_value(8), "Number of transmission buffers (FPGA->GPU) to keep around. Note that size*count buffers are kept resident throughout the program's runtime.")
    ("buffersGPU", value<unsigned>(&nbuffersGPU)->default_value(8), "Number of transmission buffers (GPU->ResultProcessing) to keep around. The required buffer size is calculated by the program and depends on the selected method and the FPGA buffer size.")
    ("timeout", value<size_t>(&timeout)->default_value(30000), "Timeout for FPGA transmissions (in ms)")
    ;
#else // GPU only
    opts_hidden_hybrid.add_options()
    ("buffersize", value<size_t>(&buffersize)->default_value(256*1024*1024), "Size for transmission buffers (GPU->ResultProcessing) in bytes.")
    ("buffersGPU", value<unsigned>(&nbuffersGPU)->default_value(16), "Number of transmission buffers (GPU->ResultProcessing) to keep around.")
    ;
#endif
#endif

    opts_positional.add("snpfile", 1);

    parse(argc, argv);
}

void Args::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular).add(opts_hidden);

    all_options.add(opts_regular_hybrid).add(opts_hidden_hybrid);

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).positional(opts_positional).run(), vars);
    notify(vars);
}

void Args::parseVars() {

    // set bools
	if (vars.count("ld") || ldfilter > 0.0)
		ld = true;
    if (vars.count("decompose"))
        decompose = true;
    if (vars.count("debug"))
        debug = true;
    if (vars.count("yaml"))
        yaml = true;
    if (vars.count("snpindex"))
        snpindex = true;

    if (!statfile.empty()) {
        bool ok = true;
        size_t pos = statfile.rfind(".bed");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        pos = statfile.rfind(".bim");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        pos = statfile.rfind(".fam");
        if (pos != string::npos && pos == statfile.size()-4)
            ok = false;
        if(!ok) {
            // Due to a potential overwrite of significant data, this error message is not written in the status file, but printed only.
            cerr << "ERROR: Status file must not end with .bed/.bim/.fam" << endl;
            exit(EXIT_FAILURE);
        }
        StatusFile::updateFile(statfile, yaml);
    }

    // set bed/bim/fam and base filename consistently
    if (snpfilebase.empty()) { // user provided --bed/bim/fam
        // check if files have the correct ending
        bool ok = true;
        size_t pos = fam.rfind(".fam");
        if (pos == string::npos || pos != fam.size()-4)
            ok = false;
        pos = bim.rfind(".bim");
        if (pos == string::npos || pos != bim.size()-4)
            ok = false;
        pos = bed.rfind(".bed");
        if (pos == string::npos || pos != bed.size()-4)
            ok = false;
        // assign snpfilebase according to provided bed file
        if (ok)
            snpfilebase.assign(bed.substr(0,pos));
        else {
            StatusFile::addError(".bed/.bim/.fam files need to have the correct endings!");
            exit(EXIT_FAILURE);
        }
    } else { // user provided filename base
        // set bed/bim/fam accordingly
        bed.assign(snpfilebase+".bed");
        bim.assign(snpfilebase+".bim");
        fam.assign(snpfilebase+".fam");
    }

    if (outprefix.empty()) { // automatically assign output prefix
        outprefix.assign(snpfilebase);
    }

    // check usage of excludeset
    if (exset_str.compare("none") && exsetA_str.compare("none")) {
        StatusFile::addError("--excludeset cannot be used together with --excludesetA");
        exit(EXIT_FAILURE);
    }
    if (exsetA_str.compare("none")) { // using exclude set only for A
        exset_str = exsetA_str;
        useexsetA = true;
    }

    // check if threshold was applied without setting nresults -> cap at 1 million best results
    if (vars.count("threshold") && !vars.count("nresults"))
    	nresults = 1000000ul;

}


bool Args::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    found |= !!this->opts_regular_hybrid.find_nothrow(optname, false);
    found |= !!this->opts_hidden_hybrid.find_nothrow(optname, false);
    return found;
}

void Args::printHelp(const string &progname, ostream &out) const {
    out << "Usage: " << progname << " <SNP input file> [options]" << endl << endl;
    out << opts_regular << endl;
#ifdef USE_CUDA_GPU
    out << opts_regular_hybrid << endl;
#endif
#ifndef NDEBUG
    out << opts_hidden << endl;
#ifdef USE_CUDA_GPU
    out << opts_hidden_hybrid << endl;
#endif
#endif
    out << endl;

    out << "The following methods are available for the -m/--method option:" << endl;
    for(const auto& m: Method::getMethods()) {
        if(m.type != Method::Type::INVALID
        		&& m.type != Method::Type::LD && m.type != Method::Type::LD3) // LD is handled seperately with the --ld switch but internally handled as an additional method
            out << "  " << left << setw(12) << m.shortName << "" << m.descriptiveName << endl;
    }
    out << endl;

    printVersion();
}

/* static */
void Args::printVersion() {
    cout << "This is version " << prog_version << ", compiled on " << prog_timestamp << endl;
    cout << "Send bugs to " << prog_bugaddress << endl;
}
