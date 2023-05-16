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

#ifndef DEVICECATEGORY_H
#define DEVICECATEGORY_H

#include <system_error>

namespace hybridsys {

using namespace std;

/**
 * @brief An error category to encapsulate runtime errors from the native Alpha Data FPGA API
 */
class FPGACategory : public error_category
{
public:
    FPGACategory() : error_category() {}

    // error_category interface
public:
    const char *name() const noexcept override ;
    string message(int code) const override;
};

/**
 * @brief An error category to encapsulate runtime errors from the native CUDA API
 */
class GPUCategory : public error_category
{
public:
    GPUCategory() : error_category() {}

    // error_category interface
public:
    const char *name() const noexcept override;
    string message(int code) const override;
};

}
#endif // DEVICECATEGORY_H
