/*
 *  ITCFileIO.h
 *  Acquirino
 *
 *  Created by CJ Bell on 5/10/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef INCLUDED_ITCFileIO_h
#define INCLUDED_ITCFileIO_h

#include "ITCAcq.h"
#include "arch_types.h"


#ifdef __cplusplus

bool WriteStructure(FILE* file, const ITCFileHeader* data);
bool WriteStructure(FILE* file, const ITCEpochHeader* data);
bool WriteStructure(FILE* file, const ITCSummedParameters* data);
bool WriteStructure(FILE* file, const ITCSineParameters* data);
bool WriteStructure(FILE* file, const ITCNoiseParameters* data);
bool WriteStructure(FILE* file, const ITCBinaryNoiseParameters* data);
bool WriteStructure(FILE* file, const ITCSegmentedNoiseParameters* data);
bool WriteStructure(FILE* file, const ITCPulseParameters* data);
bool WriteStructure(FILE* file, const ITCRampParameters* data);
bool WriteStructure(FILE* file, const ITCGenericParameters* data);
bool WriteStructure(FILE* file, const ITCGClampParameters* data);
bool WriteStructure(FILE* file, const ITCCirclePulseParameters* data);


bool ReadStructure(FILE* file, ITCFileHeader* data);
bool ReadStructure(FILE* file, ITCEpochHeader* data, int fileVersion);
bool ReadStructure(FILE* file, ITCSummedParameters* data);
bool ReadStructure(FILE* file, ITCSineParameters* data);
bool ReadStructure(FILE* file, ITCNoiseParameters* data);
bool ReadStructure(FILE* file, ITCBinaryNoiseParameters* data);
bool ReadStructure(FILE* file, ITCSegmentedNoiseParameters* data);
bool ReadStructure(FILE* file, ITCPulseParameters* data);
bool ReadStructure(FILE* file, ITCRampParameters* data);
bool ReadStructure(FILE* file, ITCGenericParameters* data);
bool ReadStructure(FILE* file, ITCGClampParameters* data);
bool ReadStructure(FILE* file, ITCCirclePulseParameters* data);

extern "C" {
#else
#include <stdbool.h>
#endif


bool WriteArray_float(FILE* file, const float* data, size_t count);
bool WriteArray_long(FILE* file, const long* data, size_t count);
bool WriteArray_int(FILE* file, const int* data, size_t count);
bool WriteArray_short(FILE* file, const short* data, size_t count);
bool Write_ITCFileHeader(FILE* file, const ITCFileHeader* data);
bool Write_ITCEpochHeader(FILE* file, const ITCEpochHeader* data);
bool Write_ITCSummedParameters(FILE* file, const ITCSummedParameters* data);
bool Write_ITCSineParameters(FILE* file, const ITCSineParameters* data);
bool Write_ITCNoiseParameters(FILE* file, const ITCNoiseParameters* data);
bool Write_ITCBinaryNoiseParameters(FILE* file, const ITCBinaryNoiseParameters* data);
bool Write_ITCSegmentedNoiseParameters(FILE* file, const ITCSegmentedNoiseParameters* data);
bool Write_ITCPulseParameters(FILE* file, const ITCPulseParameters* data);
bool Write_ITCRampParameters(FILE* file, const ITCRampParameters* data);
bool Write_ITCGenericParameters(FILE* file, const ITCGenericParameters* data);
bool Write_ITCGClampParameters(FILE* file, const ITCGClampParameters* data);
bool Write_ITCCirclePulseParameters(FILE* file, const ITCCirclePulseParameters* data);
bool Write_ITCStimulusParameters(FILE* file, const void* data, char stimulusType);


//******************************************************************************
// READING
//******************************************************************************

bool ReadArray_float(FILE* file, float* data, size_t count);
bool ReadArray_long(FILE* file, long* data, size_t count);
bool ReadArray_int(FILE* file, int* data, size_t count);
bool ReadArray_short(FILE* file, short* data, size_t count);
bool Read_ITCFileHeader(FILE* file, ITCFileHeader* data);
bool Read_ITCEpochHeader(FILE* file, ITCEpochHeader* data, int fileVersion);
bool Read_ITCSummedParameters(FILE* file, ITCSummedParameters* data);
bool Read_ITCSineParameters(FILE* file, ITCSineParameters* data);
bool Read_ITCNoiseParameters(FILE* file, ITCNoiseParameters* data);
bool Read_ITCBinaryNoiseParameters(FILE* file, ITCBinaryNoiseParameters* data);
bool Read_ITCSegmentedNoiseParameters(FILE* file, ITCSegmentedNoiseParameters* data);
bool Read_ITCPulseParameters(FILE* file, ITCPulseParameters* data);
bool Read_ITCRampParameters(FILE* file, ITCRampParameters* data);
bool Read_ITCGenericParameters(FILE* file, ITCGenericParameters* data);
bool Read_ITCGClampParameters(FILE* file, ITCGClampParameters* data);
bool Read_ITCCirclePulseParameters(FILE* file, ITCCirclePulseParameters* data);
bool Read_ITCStimulusParameters(FILE* file, void* data, char stimulusType);


#ifdef __cplusplus
} // extern "C"
#endif


#endif // INCLUDED_ITCFileIO_h
