/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt and Jan Christian Kässens,
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

#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <math.h>

using namespace std;

template<typename T>
string si_binary(T num, unsigned max_int_bits = 1) {
    static const int orders = 9;
    static const char* const symbols[orders] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"
    };

    T n = num;
    int order;
    for(order = 0; order < orders; order++) {
        if(n <= max_int_bits * 1 * 1024) {
            break;
        }
        n >>= 10;
    }

    if(order < 0)
        order = 0;
    if(order >= orders)
        order = 8;

    num >>= (10*order);
    stringstream ss;
    ss << num << " " << symbols[order];
    return ss.str();
}

template<typename T>
inline T roundToMultiple(T n, T mult) {
    if (n % mult)
        return n + mult - (n % mult);
    else
        return n;
}

template<typename T>
inline T reduceToMultiple(T n, T mult) {
    return n - (n % mult);
}

template<typename T>
inline T divideRounded(T a, T b) {
    return (a + b - 1) / b;
}

#endif // UTILS_H
