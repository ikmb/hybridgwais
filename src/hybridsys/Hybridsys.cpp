/*
 *    Copyright (C) 2018-2023 by Lars Wienbrandt, Jan Christian Kässens and Luca Schwarz,
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

#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <iostream>

extern "C" {
#include <sys/types.h>
#include <dirent.h>
}

#include "BufferFactory.h"

#include "Hybridsys.h"

using namespace std;

namespace hybridsys {

struct pci_ident {
    int vendor;
    int device;
};

struct pci_device {
    int bus;
    int slot;
};

static vector<struct pci_ident> pci_ident_fpga = {
    {0x4144, 0xADB3}, // Alpha Data; ADM-PCIE-8K5
};
static const int pci_ident_gpu_vendor = 0x10DE;

#if defined(USE_CUDA_GPU) || defined(USE_AD_FPGA)

static bool operator==(const struct pci_ident &lhs, const struct pci_ident &rhs) {
    return (lhs.vendor == rhs.vendor) && (lhs.device == rhs.device);
}

static int readHex(const string& file) {
    int result;
    ifstream f(file);
    f >> hex >> result;
    f.close();
    return result;
}

static vector<struct pci_device> findDevices(const vector<struct pci_ident>& idents) {
    vector<struct pci_device> devices;

    DIR* dir = opendir("/sys/bus/pci/devices");
    struct dirent *entry = nullptr;

    while((entry = readdir(dir))) {
        if(entry->d_type == DT_LNK) {
            struct pci_ident this_device;
            stringstream ss;
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/vendor";

            this_device.vendor = readHex(ss.str());
            ss = stringstream();
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/device";
            this_device.device = readHex(ss.str());

            if(find(begin(idents), end(idents), this_device) != end(idents)) {
                int bus = 0, slot = 0;
                string token;
                istringstream name(entry->d_name);

                getline(name, token, ':'); // extract domain and skip it
                name >> hex >> bus;
                getline(name, token, ':'); // extract bus remainder and skip it
                name >> hex >> slot;


                devices.push_back( {bus, slot} );
            }
        }
    }
    return devices;
}

static vector<struct pci_device> findDevicesFromVendor(const int ident_vendor) {
    vector<struct pci_device> devices;

    DIR* dir = opendir("/sys/bus/pci/devices");
    struct dirent *entry = NULL;

    while((entry = readdir(dir))) {
        if(entry->d_type == DT_LNK) {
            int this_device_vendor;
            stringstream ss;
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/vendor";

            this_device_vendor = readHex(ss.str());

            if(this_device_vendor == ident_vendor) {
                int bus = 0, slot = 0;
                string token;
                istringstream name(entry->d_name);

                getline(name, token, ':'); // extract domain and skip it
                name >> hex >> bus;
                getline(name, token, ':'); // extract bus remainder and skip it
                name >> hex >> slot;


                devices.push_back( {bus, slot} );
            }
        }
    }
    return devices;
}
#endif

Hybridsys::Hybridsys(vector<unsigned> allowed_fpgas __attribute__ ((unused)), vector<unsigned> allowed_gpus __attribute__ ((unused)))
    : fpgas(),
      gpus()
{
    vector<struct pci_device> devices;

    // Alpha Data API does not need to be initialized
#ifdef USE_AD_FPGA
    devices = findDevices(pci_ident_fpga);

    for(const pci_device& d: devices) {
        if(find(begin(allowed_fpgas), end(allowed_fpgas), FPGA::findIndex(d.bus, d.slot)) != end(allowed_fpgas))
            fpgas.emplace_back(d.bus, d.slot);
    }
#endif

#ifdef USE_CUDA_GPU
    // Initialize NVML
    nvmlReturn_t ret = nvmlInit();
    if (ret != NVML_SUCCESS)
        throw runtime_error(nvmlErrorString(ret));
    devices.clear();
    devices = findDevicesFromVendor(pci_ident_gpu_vendor);

    // Filter GPUs
    for(const pci_device& d: devices) {
    	gpus.emplace_back(d.bus, d.slot);
    	if(find(begin(allowed_gpus), end(allowed_gpus), gpus.back().getIndex()) == end(allowed_gpus))
    		// element is not found in the allowed list -> erase
    		gpus.pop_back();
    }
    // some devices register more than once, this ensures only one instance
    sort(gpus.begin(), gpus.end() );
	gpus.erase(unique(gpus.begin(), gpus.end()), gpus.end());
#endif
}

}
