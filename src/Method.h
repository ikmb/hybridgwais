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

#ifndef METHOD_H
#define METHOD_H

#include <istream>
#include <ostream>
#include <string>
#include <array>
#include <utility>

using namespace std;

class Method {
public:
    enum Type {
        IG,
        IG3,
        MI,
        MI3,
        BOOST,
        LOGLIN,
        LOGREG,
        LOGREG3,
        LD,
        LD3,
        INVALID,
        TYPE_MAX
    };

    struct MethodDescription {
        Method::Type type;              ///< Type ID, see Method::Type
        unsigned order;                 ///< Order of interaction (also: number of ID fields in result set)
        unsigned numScoreFields;        ///< Number of score fields in result set
        unsigned numAddScoreFields;     ///< Number of *additional* score fields when decomposed results are requested
        const char* shortName;          ///< A short identifier used for command-line args parsing and file extensions
        const char* descriptiveName;    ///< A description of the method in use, for help output
        const char** scoreFieldNames;   ///< Names (identifier) of each score field
        const char** addScoreFieldNames;///< Additional score field names for decomposed results.
    };

    static const int maxScoreComponents = 100; // must be greater than the real maximum number of components

    static const array<MethodDescription, TYPE_MAX>& getMethods();

    Method() : type(INVALID), details(false) {}
    explicit Method(Type t) : type(t), details(false) {}
    unsigned getOrder() const;
    unsigned getScoreFields() const;
    static unsigned getScoreFields(const vector<Method> &methods);
    const char* getShortName() const;
    const char* getDescription() const;
    vector<const char*> getScoreFieldNames() const;
    const char* getScoreFieldName(unsigned scorefield) const;
    Type getType() const { return type; }

    void setDetails(bool details_) { details = details_; }
    bool isDetailedComputation() const { return details; }

    static void printHeader(const vector<Method> &methods, bool snpindex, ostream &out);

private:
    Type type;
    bool details;
};

istream &operator>>(istream &in, Method &method);
ostream &operator<<(ostream &out, const Method &method);
bool operator==(const Method& lhs, const Method& rhs);
bool operator==(const Method& lhs, Method::Type t);


#endif // METHOD_H
