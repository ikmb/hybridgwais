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

inline int64_t increaseToMultipleInt(int64_t n, int64_t mult) {
	if (n < 0)
		return -(-n - ((-n) % mult));
	if (n % mult)
		return n + abs(mult) - (n % mult);
	else
		return n;
}

inline double increaseToMultipleDouble(double n, double mult) {
	if (n < 0)
		return -(-n - fmod(-n, mult));
	if (fmod(n, mult) > 0 || fmod(n, mult) < 0)
		return n + abs(mult) - fmod(n, mult);
	else
		return n;
}

inline int64_t reduceToMultipleInt(int64_t n, int64_t mult) {
	if (n < 0)
		return -(-n + abs(mult) - ((-n) % mult));
	else
		return n - (n % mult);
}

inline double reduceToMultipleDouble(double n, double mult) {
	if (n < 0)
		return -(-n + abs(mult) - fmod(-n, mult));
	else
		return n - fmod(n, mult);
}


inline int64_t divideRounded(int64_t a, int64_t b) {
	if (a < 0)
		return -(-a + abs(b) - 1) / b;
	else
		return (a + abs(b) - 1) / b;
}

#endif // UTILS_H
