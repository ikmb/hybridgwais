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
#include <algorithm>
#include "Method.h"

using namespace std;

static const char* SF_IG[] = {"IG"};
static const char* ASF_IG[] = {"MI_AB", "MI_A", "MI_B"};

static const char* SF_IG3[] = {"IG"};
static const char* ASF_IG3[] = {"MI_ABC", "MI_AB", "MI_AC", "MI_BC", "MI_A", "MI_B", "MI_C"};

static const char* SF_MI[] = {"MI"};
static const char* ASF_MI[] = {"H_ABY", "H_AB"};

static const char* SF_MI3[] = {"MI"};
static const char* ASF_MI3[] = {"H_ABCY", "H_ABC"};

static const char* SF_BOOST[] = {"BOOST_CHISQ", "BOOST_ERR", "BOOST_P-VAL"};
static const char* ASF_BOOST[] = {"BOOST_NITER", "BOOST_KSA"};

static const char* SF_LOGLIN[] = {"LL_CHISQ", "LL_ERR", "LL_P-VAL"};
static const char* ASF_LOGLIN[] = {"LL_NITER"};

static const char* SF_LOGREG[] = {"LOGR_CHISQ", "LOGR_OR", "LOGR_P-VAL", "LOGR_BETA", "LOGR_EPS"};
static const char* ASF_LOGREG[] = {"LOGR_NITER", "LOGR_DELTA", "LOGR_MINDELTA", "LOGR_RETURN"};

static const char* SF_LOGREG3[] = {"LOGR_CHISQ", "LOGR_OR", "LOGR_P-VAL", "LOGR_BETA", "LOGR_EPS"};
static const char* ASF_LOGREG3[] = {"LOGR_NITER", "LOGR_DELTA", "LOGR_MINDELTA", "LOGR_RETURN"};

static const char* SF_LD[] = {"R2_H", "R2_L"};
static const char* SF_LD3[] = {"R2_AB_H", "R2_AC_H", "R2_BC_H"};

// 1. method type
// 2. order
// 3. score fields without decomposition
// 4. additional score fields if decomposition is used
// 5. short name
// 6. description
// 7. score field names
// 8. additional score field names
static const array<Method::MethodDescription, Method::TYPE_MAX> methods {{
    { Method::Type::IG,      2, 1, 3, "ig",      "Information Gain (2-way)",                     SF_IG     , ASF_IG     },
    { Method::Type::IG3,     3, 1, 7, "ig3",     "Information Gain (3-way)",                     SF_IG3    , ASF_IG3    },
    { Method::Type::MI,      2, 1, 2, "mi",      "Mutual Information (2-way)",                   SF_MI     , ASF_MI     },
    { Method::Type::MI3,     3, 1, 2, "mi3",     "Mutual Information (3-way)",                   SF_MI3    , ASF_MI3    },
    { Method::Type::BOOST,   2, 3, 2, "boost",   "BOOST (2-way)",                                SF_BOOST  , ASF_BOOST  },
    { Method::Type::LOGLIN,  2, 3, 1, "loglin",  "Log-linear regression approximation (2-way)",  SF_LOGLIN , ASF_LOGLIN },
    { Method::Type::LOGREG,  2, 5, 4, "logreg",  "Logistic regression (2-way)",                  SF_LOGREG , ASF_LOGREG },
    { Method::Type::LOGREG3, 3, 5, 4, "logreg3", "Logistic regression (3-way)",                  SF_LOGREG3, ASF_LOGREG3},
    { Method::Type::LD,      2, 2, 0, "ld",      "Linkage disequilibrium (2-way)",               SF_LD     , NULL       },
    { Method::Type::LD3,     3, 3, 0, "ld3",     "Linkage disequilibrium between pairs (3-way)", SF_LD3    , NULL       },
    { Method::Type::INVALID, 0, 0, 0, "invalid", "(invalid)",                                    NULL      , NULL       }
}};


static const string SNPIDHEADER_2WAY = "POS A\tPOS B\tSNPID A\tSNPID B";
static const string SNPIDHEADER_3WAY = "POS A\tPOS B\tPOS C\tSNPID A\tSNPID B\tSNPID C";
static const string SNPINDEXHEADER_2WAY = "\tIDX A\tIDX B";
static const string SNPINDEXHEADER_3WAY = "\tIDX A\tIDX B\tIDX C";

const array<Method::MethodDescription, Method::Type::TYPE_MAX>& Method::getMethods() {
    return methods;
}

unsigned Method::getOrder() const {
    return methods[type].order;
}

unsigned Method::getScoreFields() const {
    unsigned ret =  methods[type].numScoreFields + (details ? methods[type].numAddScoreFields : 0);
    return ret;
}

/*static*/ unsigned Method::getScoreFields(const vector<Method> &methods_) {
    unsigned nfields = 0;
    for (const auto& m : methods_)
        nfields += m.getScoreFields();
    return nfields;
}

const char* Method::getShortName() const {
    return methods[type].shortName;
}

const char* Method::getDescription() const {
    return methods[type].descriptiveName;
}

vector<const char*> Method::getScoreFieldNames() const {
	vector<const char*> sfn;
	for (unsigned f = 0; f < methods[type].numScoreFields; f++) {
		sfn.push_back(methods[type].scoreFieldNames[f]);
	}
	if (details) {
		for (unsigned f = 0; f < methods[type].numAddScoreFields; f++) {
			sfn.push_back(methods[type].addScoreFieldNames[f]);
		}
	}
	return sfn;
}

const char* Method::getScoreFieldName(unsigned sf) const {
	if (sf >= methods[type].numScoreFields)
		return methods[type].addScoreFieldNames[sf-methods[type].numScoreFields];
	else
		return methods[type].scoreFieldNames[sf];
}

/*static*/ void Method::printHeader(const vector<Method> &methods_, bool snpindex, ostream &out) {
    if (methods_[0].getOrder() == 2) {
        out << SNPIDHEADER_2WAY;
        if (snpindex)
        	out << SNPINDEXHEADER_2WAY;
    } else {
        out << SNPIDHEADER_3WAY;
        if (snpindex)
        	out << SNPINDEXHEADER_3WAY;
    }
    for (const auto &m : methods_) {
    	for (unsigned sf = 0; sf < methods[m.type].numScoreFields; sf++)
    		out << "\t" << methods[m.type].scoreFieldNames[sf];
        if (m.details) {
        	for (unsigned sf = 0; sf < methods[m.type].numAddScoreFields; sf++)
				out << "\t" << methods[m.type].addScoreFieldNames[sf];
        }
    }
    out << endl;
}

// stream parser for kernel enum
istream& operator>>(istream &in, Method &method) {

    string token;
    in >> token;

    auto it = find_if(begin(methods), end(methods), [&](const Method::MethodDescription& m) { return m.shortName == token; });
    if(it == end(methods))
        in.setstate(ios_base::failbit);
    else
        method = Method(it->type);

    return in;
}

ostream& operator<<(ostream &out, const Method &method) {
    out << method.getShortName();
    return out;
}

bool operator==(const Method& lhs, const Method& rhs) { return lhs.getType() == rhs.getType(); }
bool operator==(const Method& lhs, Method::Type t) { return lhs.getType() == t; }

