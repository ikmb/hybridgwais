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

#ifndef SNPDB_H
#define SNPDB_H

#include <array>
#include <vector>
#include <utility>
#include <cstdint>
#include <string>

using namespace std;

// 0-based, 1st including, 2nd excluding
typedef pair<size_t,size_t> snprange;

class SNPDB
{

public:
    using AlleleCounterType = uint16_t;

    // ATTENTION!!!
    // If you change the encoding, you need to change the Boolean condition in addSampleInt
    // to add a new genotype to the buffer correctly!
    enum Genotype {
        HomozygousWild = 0,
        HomozygousVariant = 2,
        Heterozygous = 1,
        Invalid = 3
    };
    const Genotype buffer_init_gt = Invalid;
    const unsigned char buffer_init_value =
            (buffer_init_gt << 6) |
            (buffer_init_gt << 4) |
            (buffer_init_gt << 2) |
             buffer_init_gt;

    enum SexCode {
        Male,
        Female,
        Unknown
    };

    enum Phenotype {
        Control,
        Case,
        Missing
    };

    struct SampleInfo {
        string family;
        string within_family;
        string father;
        string mother;
        SexCode sex;
        Phenotype pheno;
    };

    struct SNPInfo {
        string chromosome;
        string variant_id;
        double pos_cm;
        size_t pos_bp;
        array<string, 2> alleles;
    };

    struct RegionInfo {
        // user defined global ranges
        snprange setA_g;
        snprange setB_g;
        snprange setC_g;
        snprange exset_g;
        bool useexset = false;  // user has defined an exclude set? (may be set although exset is empty!)
        bool useexsetA = false; // user defined exclude set is only for A?
        size_t excluderange; // user defined a neighboring range of SNPs to be excluded (in bp)

        // local ranges divided for processing scheme
        // the following is yet tailored for 2way interactions
        snprange l0; // part (either of A or B) left of first overlap of A and B
        snprange l1; // first overlap of A and B
        snprange l2; // part of B right of first overlap from A and B and left of second overlap from A and B
        snprange l3; // second overlap of A and B
        snprange l4; // part (either of A or B) right of second overlap of A and B
        bool l0isA = true; // indicates whether l0 is part of A (true) or part of B (false)
        bool l4isA = true; // indicates whether l4 is part of A (true) or part of B (false)
        // for 3way interactions we risk doing double calculations but we probably need these two additional intervals
        snprange l0b; // if l0 belongs to A and a possible gap in A is only covered by C (not B), l0 is the part left of the gap and l0b is right of the gap
        snprange l4b; // if l4 belongs to A and a possible gap in A is only covered by C (not B), l4 is the part left of the gap and l4b is right of the gap
        snprange lc;  // this is the local setC (if an exclude set was provided this still is a single range in the local mapping)
    };

    struct SHMHeader {
        uint64_t genotype_size;
        uint64_t case_count;
        uint64_t control_count;
        uint64_t snp_count;
    };

    template<typename T>
    class AlleleInfo {
    public:

        AlleleInfo() : counters{0,0,0,0,0,0} {}
        explicit AlleleInfo(T *counters_) {
            memcpy(this->counters, counters_, sizeof(T)*6);
        }

        inline void inc(bool is_case, SNPDB::Genotype gt, T val = 1) {
        	if (gt != SNPDB::Genotype::Invalid)
        		counters[gt + (is_case? 3 : 0)] += val;
        }

        void import(T *source) {
            memcpy(source, counters, 6 * sizeof(T));
        }

        const T *get() const {
            return counters;
        }


    private:
        T counters[6];
    };

    SNPDB(unsigned minCases = 0, unsigned minControls = 0, bool debug = false);

    SNPDB(const SNPDB&) = delete; // copy c'tor
    SNPDB& operator=(const SNPDB&) = delete; // copy assignment

    SNPDB(SNPDB&&) = default; // move c'tor
    SNPDB& operator=(SNPDB&&) & = default; // move assignment

    ~SNPDB();

    // usersetC has to be "none" for 2way tests!
    void setUserRegions(const string &usersetA_, const string &usersetB_, const string &usersetC_, const string &userexset_, bool useexsetA_, size_t excluderange_, unsigned order_) {
        usersetA  = usersetA_;
        usersetB  = usersetB_;
        usersetC  = usersetC_;
        userexset = userexset_;
        region.useexset = userexset.compare("none");
        region.useexsetA = useexsetA_;
        region.excluderange = excluderange_;
        order = order_;
        // NOTE: parsing the strings to RegionInfo is done during initialize()
    }

    void setSNPInfo(vector<SNPInfo> &&snpinfo) {
        snp_info = snpinfo;
        num_snps_global = snp_info.size();
    }

    void setSampleInfo(vector<SampleInfo> &&sampleinfo_, size_t num_cases_) {
        sample_info = sampleinfo_;
        num_cases = num_cases_;
        num_controls = sample_info.size() - num_cases_;
    }

    // reserve memory for SNP data and apply regions
    void initialize(size_t word_size);

    unsigned char *data() { return buffer; }
    const unsigned char *data() const { return buffer; }

    size_t getCaseCount() const { return num_cases; }
    size_t getControlCount() const { return num_controls; }
    size_t getSampleCount() const { return num_cases + num_controls; }
    size_t getCaseCountPadded() const { return num_cases_padded; }
    size_t getControlCountPadded() const { return num_controls_padded; }
    size_t getSampleCountPadded() const { return num_cases_padded + num_controls_padded; }
    size_t getGlobalSNPCount() const { return num_snps_global; }
    size_t getLocalSNPCount() const { return num_snps_local; }

    const vector<SampleInfo> &getSampleInfos() const { return sample_info; }
    const SampleInfo &getSampleInfo(unsigned local_index) const { return sample_info[local_index]; }
    const SNPInfo &getSNPInfo(unsigned global_index) const { return snp_info[global_index]; }
    const vector<bool> &getIgnoredSNPs() const { return ignored_snps; }
    const RegionInfo &getRegionInfo() const { return region; }
    const vector<size_t> &getMappingLtoG() const { return mapLtoG; }

    size_t getSNPSize() const { return snp_size; }
    size_t getBufferSize() const { return buffer_size; }

    inline void addSample(SNPDB::Genotype genotype, bool is_case) {
        //current_snp_alleles.inc(is_case, genotype);
        addSampleInt(genotype, (is_case ? current_case_offset : current_ctrl_offset));
    }

    void finishSNP();

    // Parses and translates a string of the form "chrN:bpstart-bpend".
    // Alternatively, a range string of the form "x-y" is allowed, where x and y are 1-based SNP indices, both including!
    // Single elements and open ranges of the form "x-" are allowed as well as the strings "all" and "none".
    // The returned pair is a pair of indices, developer-friendly 0-based and forms an interval where its left limit is included
    // and its right limit is excluded.
    snprange parseSetRange(const string &rset);

    const unsigned char * operator[](size_t snp_index) const {
        return buffer + (snp_size * snp_index);
    }

    static bool isEmpty(const snprange &r);
    static size_t lengthOf(const snprange &r);
    static bool isLeftOf(const snprange &r1, const snprange &r2); // NOTE: empty ranges are always NOT left of another range (also both ranges empty -> false)
    static bool isIncluded(const snprange &r, size_t i);
    static snprange intersec(const snprange &r1, const snprange &r2);
    static pair<snprange,snprange> cut(const snprange &r1, const snprange &r2);
    static pair<snprange,snprange> unite(const snprange &r1, const snprange &r2);
    static vector<snprange> unite(const vector<snprange> &ranges);
    static void sort(vector<snprange> &ranges);
    static vector<snprange> extract(const vector<snprange> &sets, size_t start, size_t length); // retrurns a vector of ranges that reflects the extraction of 'length' SNPs starting from SNP 'start' as if all SNPs in 'sets' were continously lined up

    static size_t getPairsCnt(const snprange &r1, const snprange &r2);
    static size_t getTriplesCnt(const snprange &r1, const snprange &r2, const snprange &r3);

    // TODO use operator<<
    static string printSNPRange(const snprange &r);

private:

    static size_t getPairsCnt(size_t a0, size_t a1, size_t b0, size_t b1);

    void applyRegions();

    // tries to find the SNP provided as chr:bp
    // returns true, if found exactly
    // returns false, if not found
    // startidx is the index in the SNP vector where the search should start
    // idx is the 0-based index value of the exact found or the next SNP if not found
    // (this may be num_snps e.g. if the chromosome is not found at all)
    bool findSNP(const string &chr, size_t bp, size_t startidx, size_t &idx);

    vector<SNPInfo> snp_info;
    vector<SampleInfo> sample_info;

    // regions
    string usersetA;
    string usersetB;
    string usersetC;
    string userexset;
    unsigned order;
    RegionInfo region;
    vector<bool> ignored_snps;
    vector<size_t> mapLtoG;

    size_t num_snps_global;
    size_t num_snps_local;
    size_t num_cases;
    size_t num_controls;
    size_t num_cases_padded;
    size_t num_controls_padded;
    size_t word_size;
    size_t snp_size;
    size_t case_size;
    typedef struct { unsigned char *bufptr; size_t bit;} genotype_offset;
    genotype_offset current_case_offset;
    genotype_offset current_ctrl_offset;
    size_t current_init_snp;
    //AlleleInfo<AlleleCounterType> current_snp_alleles;
    bool allow_buffer_resize;
    size_t buffer_size;

    const size_t snp_batch_size = 1024; // allocate this number of SNPs at once
    unsigned char *buffer;
    int buffer_fd = -1;
    size_t min_cases;
    size_t min_controls;

    bool debug;

    inline void addSampleInt(SNPDB::Genotype genotype, genotype_offset& current_offset) {

        // this requires each cell to be pre-initialized with 0xFF (i.e. an invalid genotype)
        // ATTENTION! If you change the encoding of an invalid genotype, this need to be changed as well!
        *(current_offset.bufptr) ^= ((3 ^ genotype) << current_offset.bit);

        // byte advance?
        current_offset.bit += 2;
        if(current_offset.bit == 8) {
            current_offset.bit = 0;
            ++(current_offset.bufptr);
        }
    }

    inline void resetCurrOffsets() {
        current_init_snp = 0;
        // set case offset
        current_case_offset.bufptr = buffer;
        current_case_offset.bit = 0;
        // set ctrl offset
        current_ctrl_offset.bufptr = buffer + case_size;
        current_ctrl_offset.bit = 0;
    }

};

#endif // SNPDB_H
