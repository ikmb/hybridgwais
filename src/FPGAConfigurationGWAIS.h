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

#ifndef FPGACONFIGURATIONGWAIS_H_
#define FPGACONFIGURATIONGWAIS_H_

#include <tuple>
#include "hybridsys/FPGA.h"

using namespace std;

class FPGAConfigurationGWAIS {
public:
	enum ApplicationType {
		APP_INVALID = 0,
		APP_2WAY = 2,
		APP_3WAY = 3
	};

	FPGAConfigurationGWAIS(){}

	FPGAConfigurationGWAIS(const FPGAConfigurationGWAIS &c);

	FPGAConfigurationGWAIS(const void *buffer) {
	    if (buffer != NULL)
	        parse(buffer);
	    // TODO else somehow calculate the maximum number of genotypes limited by the available GPU memory (reduced by the table buffer and the result buffer)
	}

	void parse(const void *buffer);

	unsigned getAppID() const {
		return appID;
	}
	tuple<int,int,int> getVersion() const {
		return make_tuple(versionMajor,versionMinor,versionRevision);
	}
	unsigned long getStreamFrequency() const {
		return streamFrequency;
	}
	unsigned getNumChains() const {
		return numChains;
	}
	unsigned getNumPEPerChain() const {
		return numPEPerChain;
	}
	unsigned getTableSize() const {
		return tableSize;
	}
	unsigned getNumTableEntries() const {
		return numTableEntries;
	}
	unsigned getTableEntrySize() const {
		return tableEntrySize;
	}
	unsigned getSnpIDSize() const {
		return snpIDSize;
	}
	unsigned getSnpWordSize() const {
		return snpWordSize;
	}
	unsigned getMinimumSamples() const {
		return minSamples;
	}
	unsigned getMaximumSamples() const {
		return maxSamples;
	}
	unsigned long getMaxGenotypeCount() const {
		return maxGenotypeCount;
	}

	bool operator==(const FPGAConfigurationGWAIS& other) const;

	static constexpr int major_required = 2;
	static constexpr int minor_required = 2;

private:
	unsigned appID = APP_INVALID;
	int versionMajor = 0;
	int versionMinor = 0;
	int versionRevision = 0;
	unsigned long streamFrequency = 1; // don't use 0 to prevent division by zero
	unsigned numChains = 1;
	unsigned numPEPerChain = 1;
	unsigned tableSize = 0;
	unsigned numTableEntries = 0;
	unsigned tableEntrySize = 0;
	unsigned snpIDSize = 0;
	unsigned snpWordSize = 1; // don't use 0 to prevent modulo by zero
	unsigned minSamples = 0;
	unsigned maxSamples = 0; // infinity
	unsigned long maxGenotypeCount = 0; // infinity
};

ostream& operator<<(ostream& out, const FPGAConfigurationGWAIS &c);

#endif /* FPGACONFIGURATIONGWAIS_H_ */
