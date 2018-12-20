/*
ITCACQ.H

Header file for Instrutech ITC16 acquisition library.  

Contents:
    (1) definition for file headers
    (2) definition for data acquisition epoch header
    (3) definition for stimulus parameter structure
    (4) constant definitions for data types
    (5) function prototypes
    
Created:
    5/30/97 FMR

Revision History:
    6/3/96    FMR
        (1) added function prototypes
        (2) added definition for stimulus structure

    8/14/97 FMR
        (1) changed stimulus definitions so each stimulus
            type has its own structure with its parameters
        (2) added ITCFileFooter with list of epoch positions in file
        
    9/30/97 DLR
        (1) added new variables to ITCEpochHeader
        (2) updated to version 2
        (3) fixed bug in ITCDAC2VOLTS and ITCADC2VOLTS #defines
        (4) added Offset variable to ITCFileHeader
        
    10/29/97 DLR
        (1) updated to version 3
        (2) modified epoch header to contain output conversion factors
        (3) now compiling with 68K structure alignment 

    
    12/11/97 FMR
        (1) added prototypes for analysis routines
        (2) added AmpMode variable to EpochHeader (e.g. current clamp, voltage clamp, ...)
            NOTE: did not change header length, but used spare slots at end of header
        (3) changed ITCReadDig to be consistent with ITCReadADC (returns value read rather
            than setting value of passed variable).
        (4) added ITCWriteAmpMode prototype
        (5) added prototypes for ITCGetRandomSeed and ITCNewRandomSeed
        (6) added smoothpts parameter to ITCNoiseParams structure to allow noise
            stimuli to be smoothed (with sliding boxcar average)
        (7) added a structure for generic stimuli -- to be used for storing
            parameters of stimuli generated elsewhere
        (8) added version 1 structures for conversion program
        (9) added flag for random number generator used to file header

    6/16/98 FMR
        (1) added constant to define target enviroment for compile
            useful for memory allocation purposes, e.g. Igor and MATLAB
            
    1/1/99 FMR
        (1) added SEGMENTEDNOISE stimulus type to generate stimuli like monitor (discrete updates)
        (2) added TRIGSERIES stimulus type to put out series of triggers
        (3) added prototypes for ITCSeedRandom and V4ITCNormalDeviate to update to ran2()
            random number generator from Numerical Recipes.
        (4) added BINARYNOISE stimulus type
    
    3/30/05 FMR (Version 6)
        (1) Added GCLAMP stimulus type for conductance clamp
        (2) Added CIRCLEPULSE stimulus type for monitor
        (3) Updated ITCFileHeader and ITCEpochHeader
        
	9/07 FMR/CJB (Version 8)
		update to OSX.  
		
*/

#ifndef INCLUDED_Acq_h
#define INCLUDED_Acq_h

// system includes
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
//#include <console.h>
//#include <SIOUX.h>
#include <unistd.h>
#include <fcntl.h>
//#include <Events.h>
#include <time.h>
//#include <unix.h>
//#include <InterruptDisableLib.h>
//#include <LowMem.h>

// For toolbox access
//#include <Quickdraw.h>

// ITC driver include file
#include <ITC/ITC16.h>
#include <ITC/ITC18.h>
//#include <ITC/itcmm.h>
//#include <ITC/0acqerrors.h>

#include "UserTrigger.h"


/*******************************************************************************/

// version number
#define VERSIONNUM        8

// flag for interface in use (ITC16 or ITC18)
#define ITC16 0
#define ITC18 1
#define INTERFACEFLAG ITC18

// General purpose macros (from EJC)
#define Round(x) (floor((x) + 0.5))
#define Min(x,y) ((x)<(y) ? (x) : (y))
#define Max(x,y) ((x)>(y) ? (x) : (y))
#define Clip(x,min,max) Min((Max((x),(min))),(max))

// Access Mac toolbox random number generator
// Requires inclusion of Quickdraw.h to access QD global (from EJC).
#define SetMacRandomSeed(x) ((qd.randSeed) = x)
#define GetMacRandomSeed() (qd.randSeed)
#define MacRandom() Random()
#define MinMacRandom ((long)(-32767))
#define MaxMacRandom ((long)(32767))
#define NumMacRandom ((long)(1+MaxMacRandom-MinMacRandom))

// Extreme values to use when normal range is no good (from EJC)
#define MostPositiveSingleFloat        (1.0e10)
#define MostNegativeSingleFloat        (-1.0e10)

// constants for control of FIFO during continuous acquisition
#define READFIFONUM       4000		 // read from FIFO once it has acquired this many points 
#define PTSLEFTINFIFO     100        // number of points to leave in FIFO to insure it can tell empty from full 

// other interface control constant definitions
#define ITCMINOUT      7            
#define MINLEFT 100					// minimum number of points to leave in FIFO
#define MINXFER 4000				// transfer at least this many

// other constant definitions
#define MAXSTRING        512
#define OLDMAXSTRING    128
#define PI                3.141592
#define MICROSECS        1e-6
#define PATCHCLAMP        1
#define SUCTION            2
#define ITCSAMPINTERVAL    10        // default sampling interval 10 microsec
#define MAXSEGS         4
#define MAXOUTCHAN        4        // maximum number of output channels
#define MAXINCHAN        8        // maximum number of input channels
#define MAXUNITS        16        // max number of characters in units + 1  DLR 10/2/97
#define MAXDOCHAN        8        // number of digital output channels DLR 10/3/97
#define MAXDINCHAN        16        // number of digital input channels DLR 10/3/97
#define MAXSUM            256        // maximum number of stimuli summed for single output segment
#define MAXDIGOUTPUT    32768
#define MAXDACOUTPUT    32767
#define MINDACOUTPUT    -32767
#define MAXLEDS            4
#define NUMUSERVALUES    24

// interface dependent constant definitions
#define ITCDAC2VOLTS    (32768 / 10.24)        // DLR 9/30/97
#define ITCADC2VOLTS    (32768 / 10.24)        // DLR 9/30/97
#define ITCFIFOSIZE        16384

// random number generator flags
#define RAND            0        // rand() random number generator
#define RANDOM            1        // random() random number generator (Mac Toolbox)
#define RAN2            2        // ran2() from Numerical Recipes

// operating system flags
#define WINDOWS            0
#define MACOS            1
#define UNIX            2

// stimulus types
#define SINUSOID        1
#define GAUSSIANNOISE    2
#define PULSE            3
#define RAMP            4
#define FILEWAVE        5    // DLR 10/2/97
#define GENERIC            6    // FMR 12/19/97
#define SUM                7    // FMR 12/26/97
#define TRIGSERIES        8    // FMR 11/6/98
#define SEGMENTEDNOISE  9    // FMR 1/1/99
#define BINARYNOISE        10    // FMR 1/1/99
#define GCLAMP            11  // FMR 3/05
#define CIRCLEPULSE        12    // FMR 3/05
#define ALTERNATEVOLTAGE   13    // FMR 1/09

// data storage formats for files written
#define ITCINTERLEAVEDFORMAT    0
#define ITCUNINTERLEAVEDFORMAT    1

// ITC error codes
#define ITCINITERROR    1
#define ITCSTMERROR        2
#define ITCHEADERERROR    3
#define ITCFIFOERROR    4
#define MALLOCERROR        5
#define ITCINPUTERROR    6
#define ITCFILEERROR    7
#define    ITCVERSIONERROR    8
#define ITCVERSIONWARNING 9

#define LINEARSTMAMP 0
#define LOGSTMAMP 1

// We choose 1000 here because it is slow enough for us to perform basic IO
// tasks such as reading or writing a single point without causing the ITC to
// overflow
// Historical note: a value of 10 usec) was originally used, but this caused an
// overflow condition to occur in some environments.
#define DEFAULT_SAMPLING_INTERVAL 1000


#define LOG_ERROR                                                              \
    do {                                                                       \
        printf("Fail-point (%s:%i)\n", __FILE__, __LINE__);                    \
        return(ITCFIFOERROR);                                                  \
    } while(0)

/*******************************************************************************/
// HEADERS

// general file information
struct ITCFileHeader {
    char    VersionNum;            // program version number 
    long    checksum;              // file checksum
    char    comment[MAXSTRING];    // user defined comment block
    char    DataFormat;            // storage format for data files written
    char    DataElementSize;       // size (in bytes) for each data element for portability
    float   DACtoVolts;            // conversion between volts out and DAC counts 
    float   ADCtoVolts;            // conversion between volts in and ADC counts 
    char    inputSkips;            // points deleted from first acquisition after FIFO initialization 
    long    NumEpochs;             // size of footer in bytes
    long    EpochPositionStart;    // start location of list of epoch positions
    float   Offset;                // final voltage/current offset (DLR 9/30/97)
    char    RndNumGenerator;       // flag for random number generator used (FMR 12/19/97)
    char    OperatingSystem;       // flag for target operating system (FMR 1/1/99)
    char    StmDataFormat;         // flag for data format for stimulus amplitude (FMR 2/23/00)
    long    FrameRate;             // monitor frame rate (FMR 3/05)
    char    spare[24];             // for all the missing things 
};
typedef struct ITCFileHeader ITCFileHeader;

// specific information about each acquisition epoch
struct ITCEpochHeader {
    float   time;                           // start time for acquisition relative to initialization
    char    comment[MAXSTRING];             // user defined comment block
    float   AmpGain[MAXSEGS];               // gain of amplifier 
    float   AmpOffset[MAXSEGS];             // amplifier offset in volts (DLR 9/30/97)
    float   ConvFactor[MAXSEGS];            // conversion factor for output (units/Volt) (DLR 10/29/97)
    char    SubExtOffset[MAXSEGS];          // 1 = subtract external offset voltage/current (DLR 9/30/97)
    char    ExtTrigger;                     // external trigger flag 
    long    SamplingInterval;               // sampling interval in microsec
    long    NumDataPts;                     // number of points sampled in each segment 
    long    NumSegments;                    // number of interleaved data segments to sample
    int     instructions[MAXSEGS];          // which channels sampled 
    char    InputType[MAXSEGS];             // constant descriptor of input data in each segment sampled
    char    StmOutput;                      // flag for stimulus output
    char    OutputType[MAXSEGS];            // constant descriptor of output data for each segment
    char    OutputUnits[MAXSEGS][MAXUNITS]; // units for output segments
    char    InputUnits[MAXSEGS][MAXUNITS];  // units for input segments
    char    AmpMode[MAXSEGS];               // amplifier mode
    long    MeanOutput[MAXOUTCHAN];         // mean values of output channels (FMR 3/05)    
    float   Temperature;                    // temp reading from bath thermistor (FMR 3/05)
    short   DigitalOutput;                  // status of all digital outputs (FMR 3/05)
    short   Wavelength;                     // wavelength of stimulus (FMR 3/05)
    float   NDF[MAXLEDS];                   // NDF status (FMR 3/05)
    short   LEDSetting[MAXLEDS];            // positions of LED setting switches (FMR 3/05)
    float   UserValues[NUMUSERVALUES];      // generic values for user storage (FMR 3/05)

    // Latency: number of points between reading a point and outputing the
    // response to that point. This must be at least as great as the length
    // of the ITC's hardware-pipeline.
    // "Samples" refers to hardware samples.
    // actual-latency = gClampPointsLatency * SamplingInterval
    unsigned short  gClampPointsLatency;
    
    char    spare[22];                      // again for what's missing 
};
typedef struct ITCEpochHeader ITCEpochHeader;


/*******************************************************************************/
// OLD HEADER DEFINITIONS.  Added here only if they have changed.

// version 1 file header
struct V1ITCFileHeader {
    char    VersionNum;            // program version number 
    long    checksum;              // file checksum
    char    comment[OLDMAXSTRING];    // user defined comment block
    char    DataFormat;            // storage format for data files written
    char    DataElementSize;        // size (in bytes) for each data element for portability
    float   DACtoVolts;           // conversion between volts out and DAC counts 
    float   ADCtoVolts;           // conversion between volts in and ADC counts 
    char    inputSkips;          // points deleted from first acquisition after FIFO initialization 
    long    NumEpochs;                // size of footer in bytes
    long    EpochPositionStart;    // start location of list of epoch positions
    char    spare[40];             // for all the missing things 
};
typedef struct V1ITCFileHeader V1ITCFileHeader;

// version 5 file header
struct V5ITCFileHeader {
    char    VersionNum;            // program version number 
    long    checksum;              // file checksum
    char    comment[OLDMAXSTRING];    // user defined comment block
    char    DataFormat;            // storage format for data files written
    char    DataElementSize;        // size (in bytes) for each data element for portability
    float   DACtoVolts;           // conversion between volts out and DAC counts 
    float   ADCtoVolts;           // conversion between volts in and ADC counts 
    char    inputSkips;          // points deleted from first acquisition after FIFO initialization 
    long    NumEpochs;                // size of footer in bytes
    long    EpochPositionStart;    // start location of list of epoch positions
    float   Offset;                // final voltage/current offset (DLR 9/30/97)
    char    RndNumGenerator;        // flag for random number generator used (FMR 12/19/97)
    char    OperatingSystem;        // flag for target operating system (FMR 1/1/99)
    char    StmDataFormat;            // flag for data format for stimulus amplitude (FMR 2/23/00)
    char    spare[37];             // for all the missing things 
};
typedef struct V5ITCFileHeader V5ITCFileHeader;

// version 1 epoch header
struct V1ITCEpochHeader {
    float   time;                  // start time for acquisition relative to initialization
    char    comment[OLDMAXSTRING];    // user defined comment block
    float   AmpGain[MAXSEGS];     // gain of amplifier 
    char    ExtTrigger;            // external trigger flag 
    long    SamplingInterval;        // sampling interval in microsec
    long    NumDataPts;               // number of points sampled in each segment 
    long    NumSegments;              // number of interleaved data segments to sample
    long    instructions[MAXSEGS];  // which channels sampled 
    char    InputType[MAXSEGS];    // constant descriptor of input data in each segment sampled
    char    StmOutput;             // flag for stimulus output
    char    OutputType[MAXSEGS];    // constant descriptor of output data for each segment
    char    AmpMode[MAXSEGS];        // amplifier mode
    char    spare[20-MAXSEGS];     // again for what's missing 
};
typedef struct V1ITCEpochHeader V1ITCEpochHeader;

// version 1 noise parameters structure
struct V1ITCNoiseParameters {
    long    stdev;                  // standard deviation for noise stimuli
    long    mean;                   // mean amplitude
    long    seed;                   // seed for random number generator
    long    stmpts;                 // number of points in stimulus generated 
    long    tailpts;                // points at end to set equal to mean
    long    prepts;                 // initial points set equal to mean
    char    spare[12];
};
typedef struct V1ITCNoiseParameters V1ITCNoiseParameters;

// the version 2 epoch header
struct V2ITCEpochHeader {
    float   time;                  // start time for acquisition relative to initialization
    char    comment[OLDMAXSTRING];    // user defined comment block
    float   AmpGain[MAXSEGS];     // gain of amplifier 
    float   AmpOffset[MAXSEGS];    // amplifier offset in volts (DLR 9/30/97)
    char    SubExtOffset[MAXSEGS];    // 1 = subtract external offset voltage/current (DLR 9/30/97)
    char    ExtTrigger;            // external trigger flag 
    long    SamplingInterval;        // sampling interval in microsec
    long    NumDataPts;               // number of points sampled in each segment 
    long    NumSegments;              // number of interleaved data segments to sample
    long    instructions[MAXSEGS];  // which channels sampled 
    char    InputType[MAXSEGS];    // constant descriptor of input data in each segment sampled
    char    StmOutput;             // flag for stimulus output
    char    OutputType[MAXSEGS];    // constant descriptor of output data for each segment
    char    OutputUnits[MAXSEGS][MAXUNITS];    // units for output segments
    char    InputUnits[MAXSEGS][MAXUNITS];        // units for input segments
    char    spare[20];             // again for what's missing 
};
typedef struct V2ITCEpochHeader V2ITCEpochHeader;

// version 5 epoch header
struct V5ITCEpochHeader {
    float   time;                  // start time for acquisition relative to initialization
    char    comment[OLDMAXSTRING];    // user defined comment block
    float   AmpGain[MAXSEGS];     // gain of amplifier 
    float   AmpOffset[MAXSEGS];    // amplifier offset in volts (DLR 9/30/97)
    float   ConvFactor[MAXSEGS];    // conversion factor for output (units/Volt) (DLR 10/29/97)
    char    SubExtOffset[MAXSEGS];    // 1 = subtract external offset voltage/current (DLR 9/30/97)
    char    ExtTrigger;            // external trigger flag 
    long    SamplingInterval;        // sampling interval in microsec
    long    NumDataPts;               // number of points sampled in each segment 
    long    NumSegments;              // number of interleaved data segments to sample
    long    instructions[MAXSEGS];  // which channels sampled 
    char    InputType[MAXSEGS];    // constant descriptor of input data in each segment sampled
    char    StmOutput;             // flag for stimulus output
    char    OutputType[MAXSEGS];    // constant descriptor of output data for each segment
    char    OutputUnits[MAXSEGS][MAXUNITS];    // units for output segments
    char    InputUnits[MAXSEGS][MAXUNITS];        // units for input segments
    char    AmpMode[MAXSEGS];        // amplifier mode
    char    spare[20-MAXSEGS];     // again for what's missing 
};
typedef struct V5ITCEpochHeader V5ITCEpochHeader;

/*******************************************************************************/
// stimulus specifications

// pulse stimuli
struct ITCPulseParameters {
     long    amp;                    // pulse amp for noise stimuli
     long    mean;                   // mean amplitude
     long    stmpts;                 // number of points in stimulus generated 
     long    tailpts;                // points at end to set equal to mean
     long    prepts;                 // initial points set equal to mean
     char    spare[12];
};
typedef struct ITCPulseParameters ITCPulseParameters;

// sinusoidal stimuli
struct ITCSineParameters { 
     long    amp;                    // sine amp for noise stimuli
     long    mean;                   // mean amplitude
     long    period;                 // period in units of sampling interval
     long    phase;                  // initial phase for sine stimulus
     long    stmpts;                 // number of points in stimulus generated 
     long    tailpts;                // points at end to set equal to mean
     long    prepts;                 // initial points set equal to mean
     char    spare[12];
};
typedef struct ITCSineParameters ITCSineParameters;

// noise stimuli
struct ITCNoiseParameters {
    long    stdev;                  // standard deviation for noise stimuli
    long    mean;                   // mean amplitude
    long    seed;                   // seed for random number generator
    long    stmpts;                 // number of points in stimulus generated 
    long    tailpts;                // points at end to set equal to mean
    long    prepts;                 // initial points set equal to mean
    long    smoothpts;                 // number of points to smooth over with sliding boxcar average
    long    freqcutoff;             // frequency cutoff for smoothing (1/1/99 FMR)
    long    numfilt;                 // number of filters in cascade for smoothing (1/1/99 FMR)
};
typedef struct ITCNoiseParameters ITCNoiseParameters;

// segmented noise stimuli 1/1/99 FMR
struct ITCSegmentedNoiseParameters {
    long    stdev;                  // standard deviation for noise stimuli
    long    mean;                   // mean amplitude
    long    seed;                   // seed for random number generator
    long    stmpts;                 // number of points in stimulus generated 
    long    tailpts;                // points at end to set equal to mean
    long    prepts;                 // initial points set equal to mean
    long    segmentpts;             // number of points for each independent noise value
    char    spare[12];
};
typedef struct ITCSegmentedNoiseParameters ITCSegmentedNoiseParameters;

// binary noise stimuli 1/1/99 FMR
struct ITCBinaryNoiseParameters {
    long    amp;                       // amplitude for noise stimuli
    long    mean;                   // mean amplitude
    long    seed;                   // seed for random number generator
    long    stmpts;                 // number of points in stimulus generated 
    long    tailpts;                // points at end to set equal to mean
    long    prepts;                 // initial points set equal to mean
    long    segmentpts;             // number of points for each independent noise value
    char    spare[12];
};
typedef struct ITCBinaryNoiseParameters ITCBinaryNoiseParameters;

// ramp stimuli            FMR 12/11/97
struct ITCRampParameters {
     long    amp;                    // pulse amp for noise stimuli
     long    mean;                   // mean amplitude
     long    stmpts;                 // number of points in stimulus generated 
     long    tailpts;                // points at end to set equal to mean
     long    prepts;                 // initial points set equal to mean
     char    spare[12];
};
typedef struct ITCRampParameters ITCRampParameters;

// trigger series            FMR 11/6/98
struct ITCTrigSeriesParameters {
     long    CurDigOut;              // constant digital output to add to trigger series
     long    TrigOut;                // output corresponding to digital out to toggle (e.g. 2 = TTL Out #1)
     long    interval;               // number of points between triggers
     long    stmpts;                 // number of points in stimulus generated 
     long    tailpts;                // points at end to set equal to mean
     long    prepts;                 // initial points set equal to mean
     char    spare[12];
};
typedef struct ITCTrigSeriesParameters ITCTrigSeriesParameters;

// trigger series            FMR 1/09
struct ITCAlternateVoltageParameters {
	long    interval;               // number of points between triggers
	long	amp;					// amplitude of step
	long    mean;                   // mean amplitude
	long    stmpts;                 // number of points in stimulus generated 
	long    tailpts;                // points at end to set equal to mean
	long    prepts;                 // initial points set equal to mean	
	char    spare[12];
};
typedef struct ITCAlternateVoltageParameters ITCAlternateVoltageParameters;

// generic stimuli            FMR 12/11/97
struct ITCGenericParameters {
     long    amp;                    // pulse amp for noise stimuli
     long    mean;                   // mean amplitude
     long    stmpts;                 // number of points in stimulus generated 
     long    tailpts;                // points at end to set equal to mean
     long    prepts;                 // initial points set equal to mean
     char    spare[12];
};
typedef struct ITCGenericParameters ITCGenericParameters;

// sum of several stimuli        FMR 12/26/97
struct ITCSummedParameters {
    long    NumStimuli;            // number of stimuli summed together for this segment
    char    StmType[MAXSUM];        // constant descriptor for each stimulus type in sum
    void	*StmParams[MAXSUM];    // pointers to stimulus structures
};
typedef struct ITCSummedParameters ITCSummedParameters;

// filewave stimuli        DLR 10/2/97
struct ITCFilewaveParameters
{
    short*  data;                   // where the data resides
    long    pts;                    // number of data points
};
typedef struct ITCFilewaveParameters ITCFilewaveParameters;

// conductance clamp stimuli        FMR 3/30/05
struct ITCGClampParameters
{
    float   ExcitatoryVRev;        // reversal potential for excitatory inputs
    float   InhibitoryVRev;        // reversal potential for inhibitory inputs
};
typedef struct ITCGClampParameters ITCGClampParameters;

// pulsed circular spot stimuli from monitor    FMR 3/30/05
struct ITCCirclePulseParameters
{
    long    XCenter;               // x coordinate of center
    long    YCenter;               // y coordinate of center
    long    diameter;              // diameter of spot
    long    amp[3];                // pulse amp of each monitor gun
    long    mean[3];               // mean amplitude of each monitor gun
    long    stmframes;             // number of frames in pulse
    long    tailframes;            // frames at end to set equal to mean
    long    preframes;             // initial frames set equal to mean
};
typedef struct ITCCirclePulseParameters ITCCirclePulseParameters;

/*******************************************************************************/                                              
// function declarations

// Initializations
long ITCInit(short **InputBuffer, short **OutputBuffer, long *bufsize, int deviceNumber);
long ITCSetSeq(char *InputString, char *OutputString);
long ITCSetAcqParams(long SamplingInterval, long NumPts, char StmOutput, char ExtTrigger);
long ITCInitAcq();
long ITCFinish();
long ITCQuickReInit(int deviceNumber);

// Acquisition
long ITCAcqEpoch();
long ITCAcqEpochGClamp(double ADCCounts2mV, double DACCounts2I);
long ITCReadDig(long *ValueRead);
long ITCSetDig(short Value);
long ITCSetDAC(short Value, char Number);
long ITCReadADC(long *ValueRead, char Number);
long ITCFastAcqEpoch(short ScopeTrigger, short DigVal);

// Utilities
void ITCInterleaveData(short *InterleavedDataPtr, short **DataPtr, long SegmentCnt, long SegmentPts);
void ITCSplitData(short *InterleavedDataPtr, short **DataPtr, long SegmentCnt, long SegmentPts);
//long ITCKeyEvent();

// File access and Saving
long ITCReadFileHeader(char ElementCode, float   HeaderValue);
long ITCWriteFileHeader(char ElementCode, float   HeaderValue);
long ITCCreateFile(char *FileName, char OverWriteFlag);
long ITCSaveToFile();
void ITCCloseFile();
long ITCGetNumEpochs ();
long ITCOpenFile(char *FileName);
long ITCOpenFileForEdit(char *FileName);
long ITCGetInterval ();
void ITCWriteOffset ( float   );
void ITCWriteInputUnits ( char *units, long seg );
char *ITCGetInputUnits ( long seg );
void ITCWriteOutputUnits ( char *units, long seg );
char *ITCGetOutputUnits ( long seg );
long ITCGetDataPts (), ITCGetSegments ();
char *ITCGetEpochCom ();
char *ITCGetFileCom ();
void ITCWriteFileCom(char *comment);
void ITCWriteEpochCom(char *comment);
void ITCWriteAmpGain(float   gain, float   offset, char subtract, long SegmentNumber);
void ITCGetAmpSettings ( float   *, float   *, char *, long );
float   ITCGetOffset ();
float   ITCGetEpochTime ();
void ITCGetStimParams ( char SegmentNumber, char *Params, int size );    // DLR 10/29/97
int ITCGetStimType ( char SegmentNumber );    // DLR 10/29/97
void ITCWriteConvFactor ( float   conv_factor, long seg ); // DLR 10/29/97
float   ITCGetConvFactor ( long seg );  // DLR 10/29/97
long ITCGetSampInterval ();  // DLR 10/29/97
void ITCWriteAmpMode(long mode, long SegmentNumber); // FMR 11/29/97
long ITCRetrieveEpoch(long EpochNumber);
long ITCRetrieveEpochStm(long EpochNumber);            // FMR 12/11/97
long ITCRetrieveEpochStmGClamp(long EpochNumber, long SegmentNumber);        // FMR 4/05
long ITCEpochParameters(long EpochNumber);
long ITCRetrieveStmAmp(long *ValueRead, long StmType, void *params);        // FMR 12/11/97
long ITCRetrieveStmMean(long *ValueRead, long StmType, void *params);        // FMR 4/25/98
long ITCRetrieveStmTailPts(long *ValueRead, long StmType, void *params);    // FMR 12/11/97
long ITCRetrieveStmPts(long *ValueRead, long StmType, void *params);        // FMR 12/11/97
long ITCRetrieveStmPrePts(long *ValueRead, long StmType, void *params);        // FMR 12/11/97
long ITCRetrieveRandomSeed(long *ValueRead, long SegmentNumber);    // FMR 12/11/97
long ITCStmParamsSize(long *StmParamsSize, char StimulusType);        // FMR 12/11/97
long ITCReadStmParams(int SegmentNumber);                            // FMR 12/27/97
long ITCWriteStmParams(int SegmentNumber, long *EpochSize);            // FMR 12/27/97
void ITCWriteTemperature(float   Temp);                                // FMR 4/05
void ITCWriteMeanOutput(short OutputValue, long OutputChannel);        // FMR 4/05
void ITCWriteUserValue(float   Value, long Number);                    // FMR 4/05
void ITCWriteLEDSetting(short Setting, long LEDNumber);                // FMR 4/05
void ITCWriteNDF(float   NDFValue, long NDFNumber);                    // FMR 4/05    
void ITCWriteWavelength(short Wavelength);                            // FMR 4/05
void ITCWriteDigitalOutput(short CurrentDigitalOutput);                // FMR 4/05    

// Analysis FMR 11/29/97
long ITCInitAnalysis(short **InputBuffer, short **OutputBuffer, long *npts);
void ITCThresh(long Threshold, long NumPts, long RefractTime, float   *Buffer);
void ITCThreshMat(long Threshold, long NumPts, long RefractTime, double *Buffer);
long ITCSpikeTrigPower(long NumPts, long STAPts, long BinPoints, long StmPower, float   *stm, float   *spikes, float   *var);
long ITCSpikeTrigCovar(long NumPts, long STAPts, long BinPoints, float   *stm, float   *spikes, float   *sta, float   *covar, float   *stmcovar);
long ITCSpikeTrigAve(long NumPts, long STAPts, long BinPoints, double *stm, double *spikes, double *var);
long ITCAnalogCovariance(long NumPts, long STAPts, double *stm, double *volt, double *covar, double *stmcovar);
long ITCResponseCovariance(long NumPts, double *response, double *covar);

// Stimulus Generation
long ITCUpdateStm();
void ITCStmCreateSine(short *stm, ITCSineParameters *params);
long ITCStmCreateNoise(short *stm, ITCNoiseParameters *params);
long ITCStmCreateSegmentedNoise(short *stm, ITCSegmentedNoiseParameters *params);
void ITCStmCreatePulse(short *stm, ITCPulseParameters *params);
void ITCStmCreateRamp(short *stm, ITCRampParameters *params);
void ITCStmCreateAltVoltage(short *stm, ITCAlternateVoltageParameters *params);
void ITCStmCreateGClamp(short *stm, ITCGClampParameters *params);        // FMR 4/05
void ITCStmCreateFilewave(short *stm, ITCFilewaveParameters *params);    // DLR 10/2/97
void ITCStmCreateTrigSeries(short *stm, ITCTrigSeriesParameters *params);  // FMR 11/6/98
long ITCStmCreateBinaryNoise(short *stm, ITCBinaryNoiseParameters *params);    // FMR 1/1/99
long ITCStmCreateSum(short *stm, ITCSummedParameters *params);        // FMR 12/27/97
long ITCCalcStimulus(char StimType, short *stm, void *params);
long ITCSummedParams(char SegmentNumber, ITCSummedParameters *SummedParams);    // FMR 12/28/97
long ITCSineParams(char SegmentNumber, ITCSineParameters *SineParams);
long ITCPulseParams(char SegmentNumber, ITCPulseParameters *PulseParams);
long ITCBinaryNoiseParams(char SegmentNumber, ITCBinaryNoiseParameters *NoiseParams);    // FMR 1/1/99
long ITCRampParams(char SegmentNumber, ITCRampParameters *RampParams);
long ITCAltVoltageParams(char SegmentNumber, ITCAlternateVoltageParameters *AltVoltParams);
long ITCNoiseParams(char SegmentNumber, ITCNoiseParameters *NoiseParams);
long ITCGClampParams(char SegmentNumber, ITCGClampParameters *GClampParams, long samplesLatency);
long ITCSegmentedNoiseParams(char SegmentNumber, ITCSegmentedNoiseParameters *NoiseParams);
long ITCTrigSeriesParams(char SegmentNumber, ITCTrigSeriesParameters *TrigSeriesParams); // FMR 11/6/98
long ITCFilewaveParams(char SegmentNumber, ITCFilewaveParameters *FilewaveParams); // DLR 10/2/97
long ITCGenericParams(ITCGenericParameters *Params);                            // FMR 1/1/99
long ITCAddSineParams(char SegmentNumber, ITCSineParameters *SineParams);        // FMR 12/29/97
long ITCAddPulseParams(char SegmentNumber, ITCPulseParameters *PulseParams);    // FMR 12/29/97
long ITCAddRampParams(char SegmentNumber, ITCRampParameters *RampParams);        // FMR 12/29/97
long ITCAddNoiseParams(char SegmentNumber, ITCNoiseParameters *NoiseParams);    // FMR 12/29/97
long ITCAddSegmentedNoiseParams(char SegmentNumber, ITCSegmentedNoiseParameters *NoiseParams);    // FMR 12/29/97
long ITCAddBinaryNoiseParams(char SegmentNumber, ITCBinaryNoiseParameters *NoiseParams);    // FMR 1/1/99
long ITCAddFilewaveParams(char SegmentNumber, ITCFilewaveParameters *FilewaveParams);    // FMR 12/29/97
long ITCNewRandomSeed(long seed, long SegmentNumber, long SumNumber);
long ITCGetRandomSeed(long SegmentNumber);
void ITCFreeParams(char SegmentNumber);                                        // FMR 12/30/97
long ITCAddCirclePulseParams(char SegmentNumber, ITCCirclePulseParameters *params);        // FMR 4/05
long ITCCirclePulseParams(char SegmentNumber, ITCCirclePulseParameters *params);        // FMR 4/05
long ITCAddGClampParams(char SegmentNumber, ITCGClampParameters *GClampParams);         // FMR 6/05

// Utilities
float   V4ITCNormalDeviate();            // FMR 1/1/99
float   ITCNormalDeviate(long *idum);    // FMR 1/1/99
void WaitKey();
void ITCStartMsecTimer();            // FMR 1/14/98
long ITCReadMsecTimer();            // FMR 1/14/98
long ITCWaitMsec(long msec);        // FMR 1/14/98    
char *ITCMalloc(char PersistantFlag, long NumBytes);    // FMR 6/16/98
char ITCFree(void *MemPtr);                                            // FMR 6/16/98
void ITCSeedRandom(long *seed);        // FMR 1/1/99
void ITCfft(int type,int num,float   *dat);            // FMR 1/1/99
int ConvertMacPathToUnix(const char * inPath, char * outPath, int outPathMaxLen);

#define MIN(A,B) (((A) < (B)) ? (A) : (B))

#endif // INCLUDED_Acq_h

