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

#include <stdint.h>
#include "hybridsys/FPGA.h"
#include "FPGAConfigurationGWAIS.h"

using namespace std;

FPGAConfigurationGWAIS::FPGAConfigurationGWAIS(const FPGAConfigurationGWAIS& c) :
	appID(c.appID),
	versionMajor(c.versionMajor),
	versionMinor(c.versionMinor),
	versionRevision(c.versionRevision),
	streamFrequency(c.streamFrequency),
	numChains(c.numChains),
	numPEPerChain(c.numPEPerChain),
	tableSize(c.tableSize),
	numTableEntries(c.numTableEntries),
	tableEntrySize(c.tableEntrySize),
	snpIDSize(c.snpIDSize),
	snpWordSize(c.snpWordSize),
	minSamples(c.minSamples),
	maxSamples(c.maxSamples),
	maxGenotypeCount(c.maxGenotypeCount)
	{}

void FPGAConfigurationGWAIS::parse(const void *buffer) {
    const uint16_t *info = reinterpret_cast<const uint16_t*>(buffer);

    appID            = info[0] >> 8;
    versionMajor     = info[0] & 0xff;
    versionMinor     = info[1] >> 8;
    versionRevision  = info[1] & 0xff;
    streamFrequency  = (((uint32_t)info[3]) << 16) | ((uint32_t) info[2]);
    numChains        = info[4];
    numPEPerChain    = info[5];
    tableSize        = info[6];
    numTableEntries  = info[7];
    tableEntrySize   = info[8];
    // info[9] is reserved
    snpIDSize        = info[10];
    snpWordSize      = info[11];
    minSamples       = info[12];
    maxSamples       = info[13];
    maxGenotypeCount = 1UL << static_cast<unsigned long>(info[14]);
    // info[15] is reserved (and includes the status byte)
}

bool FPGAConfigurationGWAIS::operator==(const FPGAConfigurationGWAIS &other) const {
    return
        (appID            == other.appID           ) &&
        (versionMajor     == other.versionMajor    ) &&
        (versionMinor     == other.versionMinor    ) &&
        (versionRevision  == other.versionRevision ) &&
        (streamFrequency  == other.streamFrequency ) &&
        (numChains        == other.numChains       ) &&
        (numPEPerChain    == other.numPEPerChain   ) &&
        (tableSize        == other.tableSize       ) &&
        (numTableEntries  == other.numTableEntries ) &&
        (tableEntrySize   == other.tableEntrySize  ) &&
        (snpIDSize        == other.snpIDSize       ) &&
        (snpWordSize      == other.snpWordSize     ) &&
        (minSamples       == other.minSamples      ) &&
        (maxSamples       == other.maxSamples      ) &&
        (maxGenotypeCount == other.maxGenotypeCount);
}

ostream& operator<<(ostream& out, const FPGAConfigurationGWAIS &c) {
    out << "\tApplication type:             ";
    switch(c.getAppID()) {
    case FPGAConfigurationGWAIS::APP_2WAY: out << "2-way interactions"; break;
    case FPGAConfigurationGWAIS::APP_3WAY: out << "3-way interactions"; break;
    default: out << "unknown"; break;
    }
    out << endl;
    out << "\tVersion:                      " << get<0>(c.getVersion()) << "." << get<1>(c.getVersion()) << " rev " << get<2>(c.getVersion()) << endl;
    out << "\tFrequency of stream pipeline: " << c.getStreamFrequency() << " Hz" << endl;
    out << "\tNumber of PEs:                " << c.getNumChains() << " x " << c.getNumPEPerChain() << " = " << (c.getNumChains() * c.getNumPEPerChain()) << endl;
    out << "\tTable size:                   " << c.getTableSize() << " byte" << endl;
    out << "\tNumber of table entries:      " << c.getNumTableEntries() << endl;
    out << "\tSize of table entry:          " << c.getTableEntrySize() << " bit" << endl;
    out << "\tSize of SNP ID:               " << c.getSnpIDSize() << " bit" << endl;
    out << "\tSize of raw SNP data word:    " << c.getSnpWordSize() << " bytes" << endl;
    return out;
}
