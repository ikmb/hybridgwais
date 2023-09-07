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


// 1. method type
// 2. order
// 3. score fields without decomposition
// 4. additional score fields if decomposition is used
// 5. short name
// 6. description
static const array<Method::MethodDescription, Method::TYPE_MAX> methods {{
    { Method::Type::IG,      2, 1, 3, "ig",      "Information Gain (2-way)"},
    { Method::Type::IG3,     3, 1, 7, "ig3",     "Information Gain (3-way)"},
    { Method::Type::MI,      2, 1, 2, "mi",      "Mutual Information (2-way)"},
    { Method::Type::MI3,     3, 1, 2, "mi3",     "Mutual Information (3-way)"},
    { Method::Type::BOOST,   2, 3, 2, "boost",   "BOOST (2-way)"},
    { Method::Type::LOGLIN,  2, 3, 1, "loglin",  "Log-linear regression approximation (2-way)"},
    { Method::Type::LOGREG,  2, 5, 4, "logreg",  "Logistic regression (2-way)"},
    { Method::Type::LOGREG3, 3, 5, 4, "logreg3", "Logistic regression (3-way)"},
    { Method::Type::LD,      2, 2, 0, "ld",      "Linkage disequilibrium (2-way)"},
    { Method::Type::LD3,     3, 3, 0, "ld3",     "Linkage disequilibrium between pairs (3-way)"},
    { Method::Type::INVALID, 0, 0, 0, "invalid", "(invalid)"},
}};

static const string SNPIDHEADER_2WAY = "POS A\tPOS B\tSNPID A\tSNPID B";
static const string SNPIDHEADER_3WAY = "POS A\tPOS B\tPOS C\tSNPID A\tSNPID B\tSNPID C";
static const string SNPINDEXHEADER_2WAY = "\tIDX A\tIDX B";
static const string SNPINDEXHEADER_3WAY = "\tIDX A\tIDX B\tIDX C";

// TODO make compatible for multiple methods
static const array<string, Method::TYPE_MAX> SCOREHEADERS {
  "\tIG", // IG2
  "\tIG", // IG3
  "\tMI", // MI2
  "\tMI", // MI3
  "\tBOOST_CHISQ\tBOOST_ERR\tBOOST_P-VAL", // BOOST
  "\tLL_CHISQ\tLL_ERR\tLL_P-VAL", // LOGLIN
  "\tLOGR_CHISQ\tLOGR_OR\tLOGR_P-VAL\tLOGR_BETA\tLOGR_EPS", // LOGREG
  "\tLOGR_CHISQ\tLOGR_OR\tLOGR_P-VAL\tLOGR_BETA\tLOGR_EPS", // LOGREG3
  "\tR2_H\tR2_L", // LD2
  "\tR2_AB_H\tR2_AC_H\tR2_BC_H", // LD3
  "" // INVALID
};

// TODO make compatible for multiple methods
static const array<string, Method::TYPE_MAX> SCOREHEADER_DECOMPAPPEND {
  "\tMI_AB\tMI_A\tMI_B", // IG2
  "\tMI_ABC\tMI_AB\tMI_AC\tMI_BC\tMI_A\tMI_B\tMI_C", // IG3
  "\tH_ABY\tH_AB", // MI2
  "\tH_ABCY\tH_ABC", // MI3
  "\tBOOST_#ITER\tBOOST_KSA", // BOOST
  "\tLL_#ITER", // LOGLIN
  "\tLOGR_#ITER\tLOGR_DELTA\tLOGR_MINDELTA\tLOGR_RETURN", // LOGREG
  "\tLOGR_#ITER\tLOGR_DELTA\tLOGR_MINDELTA\tLOGR_RETURN", // LOGREG3
  "", // LD2
  "", // LD3
  "" // INVALID
};

const array<Method::MethodDescription, Method::Type::TYPE_MAX>& Method::getMethods() {
    return methods;
}

unsigned Method::getOrder() const {
    return methods[type].order;
}

unsigned Method::getScoreFields() const {
    unsigned ret =  methods[type].numScoreFields + (details ? methods[type].numScoreComponents : 0);
    return ret;
}

/*static*/ unsigned Method::getScoreFields(const vector<Method> &methods_) {
    unsigned nfields = 0;
    for (const auto& m : methods_)
        nfields += m.getScoreFields();
    return nfields;
}

const char *Method::getShortName() const {
    return methods[type].shortName;
}

const char *Method::getDescription() const {
    return methods[type].descriptiveName;
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
        out << SCOREHEADERS[m.type];
        if (m.details)
            out << SCOREHEADER_DECOMPAPPEND[m.type];
    }
    out << endl;
}

// stream parser for kernel enum
istream &operator>>(istream &in, Method &method) {

    string token;
    in >> token;

    auto it = find_if(begin(methods), end(methods), [&](const Method::MethodDescription& m) { return m.shortName == token; });
    if(it == end(methods))
        in.setstate(ios_base::failbit);
    else
        method = Method(it->type);

    return in;
}

ostream &operator<<(ostream &out, const Method &method) {
    out << method.getShortName();
    return out;
}

bool operator==(const Method& lhs, const Method& rhs) { return lhs.getType() == rhs.getType(); }
bool operator==(const Method& lhs, Method::Type t) { return lhs.getType() == t; }

