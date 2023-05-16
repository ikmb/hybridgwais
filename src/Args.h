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

#ifndef ARGS_H_
#define ARGS_H_

#include <string>
#include <utility>
#include <array>

#include <boost/program_options.hpp>

#include "Method.h"
#include "GPUEngine.h"

namespace bpo = boost::program_options;

using namespace std;

/**
 * Class for storing and retrieving command-line arguments.
 */
class Args {
public:

	static Args parseArgs(int argc, char *argv[]);

private:
    /**
     * Returns the argument with the given (long) name. The template argument
     * specifies the return type you want the parameter casted to.
     * @param name argument name
     * @return the variable value
     */
    template<typename T> T get(const string &name) const {
        auto where = vars.find(name);
        if(where == end(vars)) {
            if(!isDefined(name))
                throw invalid_argument("Option undefined: " + name + " (This is a bug)");
            else
                throw out_of_range("Option has not been specified and does not have a default value associated: " + name);
        }
        return where->second.as<T>();
    }

    /**
     * Counts the number of argument occurences. Mainly useful for boolean switches
     * and/or counting flags, such as verbosity or debug mode.
     * @param name argument name
     * @return argument value
     */
    unsigned int count(const string &name) const {
        if(!isDefined(name))
            throw invalid_argument("Option undefined: " + name + " (This is a bug)");
        return vars.count(name);
    }

    bool operator()(const string &name) const {
        return count(name) > 0;
    }
public:
    Args(Args&& other) = default;

    string outprefix;
    string statfile;
    vector<Method> methods;
    bool ld = false;
    GPUEngine::score_type ldfilter = 0.0;
    size_t nresults;
    bool decompose = false;
    GPUEngine::score_type threshold;
    string snpfilebase;
    string bed;
    string bim;
    string fam;
    string setA_str;
    string setB_str;
    string setC_str;
    string exset_str;
    bool useexsetA = false;
    size_t excluderange;

    vector<unsigned> fpgas;
    vector<unsigned> gpus;

    bool debug = false;
    bool yaml = false;

    size_t buffersize;
    unsigned nbuffersFPGA;
    unsigned nbuffersGPU;
    size_t timeout;

protected:
    /** Constructs the arguments list and adds all defined options */
    Args();
    Args(Args const &);
    void operator=(Args const &);

    void parse(int argc, char *argv[]);
    void parseVars();
    bool isDefined(const string &optname) const;

    void printHelp(const string &progname, ostream &out) const;
    static void printVersion();

    bpo::options_description opts_regular;        /**< regular options, shown on help output */
    bpo::options_description opts_regular_hybrid; /**< regular options, shown on help output */
    bpo::options_description opts_hidden;         /**< hidden options */
    bpo::options_description opts_hidden_hybrid;  /**< hidden options */
    bpo::positional_options_description opts_positional;    /**< positional options (without name) */

    bpo::variables_map vars;    /**< this is where the option values go */

    /** dump all options to the given stream */
    friend ostream &operator<<(ostream &out, const Args &args);

private:
    /** parses the main() options into the variable map */
    Args(int argc, char *argv[]);

    // used only internally
    string exsetA_str;
};

#endif /* ARGS_H_ */

