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

#ifndef DEVICE_H
#define DEVICE_H

#include <string>

namespace hybridsys {

using namespace std;

class Device {
public:
    virtual ~Device() {}

    /**
     * @brief Getter method for the logical PCI Express bus ID
     * @return the assigned bus ID
     */
    virtual int getBus() = 0;

    /**
     * @brief Getter method for the logical PCI Express slot/device ID
     * @return the assigned slot/device ID
     */
    virtual int getSlot() = 0;

    /**
     * @brief Getter method for the serial number
     * @return a string representation of the device's serial number
     * @throws runtime_error with FPGACategory and `ADMXRC3_BAD_DRIVER` error code if the device does not have a valid bitstream loaded (FPGA only)
     */
    virtual const string& getSerialNumber() = 0;

    /**
     * @brief Getter method for the device index that is used to instantiate the underlying native API
     * @return a logical device index. It is not guaranteed to be persistent among reboots.
     * @throws runtime_error with FPGACategory and `ADMXRC3_BAD_DRIVER` error code if the device does not have a valid bitstream loaded (FPGA only)
     */
    virtual int getIndex() const = 0;
};

}

#endif // DEVICE_H
