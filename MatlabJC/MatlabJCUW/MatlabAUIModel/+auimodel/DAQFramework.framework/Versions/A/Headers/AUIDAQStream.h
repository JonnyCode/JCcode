//
//  AUIDAQStreamInfo.h
//  AcqUI
//
//  Created by Barry Wark on 4/4/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIExternalDevice.h>
#import <DAQFramework/AUIIOStream.h>

#define DAQSTREAMINFO_SAMPLINGRATE_ERROR 0

//If you change this, you must also change the byte-swapping functions used in daqToAnalogData
//#define ANALOG_USER_SAMPLE_TYPE float
typedef float ANALOG_USER_SAMPLE_TYPE;


@class AUIDAQController;
/*!
    @class
    @abstract    AUIIOStream subclass for ADC/DAC streams.
    @discussion  ADC/DAC samples must be signed integers.
*/

@interface AUIDAQStream : AUIIOStream {
	NSMutableData *data;
	double			samplingRate;	// Hz
	long backgroundSample;
	
    unsigned long bytesTransfered;
    double voltToSampleMult;
    unsigned sampleBytes;
    
    StreamStatus	status;
    
    NSNumber *lastSample;
    NSMutableData *backgroundData;
}

+ (unsigned)bitsPerAnalogUserSample;


- (id)initWithIOController:(AUIDAQController*)c 
             channelNumber:(unsigned long)theChannelNumber
                      type:(StreamType)theType
               sampleBytes:(unsigned)_sampleBytes 
          voltToSampleMult:(double)_voltToSampleMult; //designated initializer

- (id)initWithIOController:(AUIIOController*)c
             channelNumber:(unsigned long)theChannelNumber
                      type:(StreamType)theType;

- (NSNumber*)background;
- (NSData*)backgroundDataForSamples:(unsigned long)nSamples;
- (void)setBackground:(NSNumber*)newBackground;

- (StreamStatus)status;
- (void)setStatus:(StreamStatus)newStatus;

- (unsigned)dataLength;

- (unsigned long)bytesTransfered;
- (void)setBytesTransfered:(unsigned long)newBytesTransfered;

- (double)voltToSampleMult;
- (void)setVoltToSampleMult:(double)newVoltToSampleMult;
- (unsigned)sampleBytes;
- (void)setSampleBytes:(unsigned)newSampleBytes;

- (NSNumber *)lastSample;
- (void)setLastSample:(NSNumber *)aLastSample;

- (NSData*)analogDataToDAQ:(NSData*)inData;
- (NSData*)analogSampleToDAQ:(ANALOG_USER_SAMPLE_TYPE)sample;
- (NSData*)daqToAnalogData:(NSData*)inData byteOrder:(CFByteOrder)byteOrder; //output is in byteOrder format
- (ANALOG_USER_SAMPLE_TYPE)daqToAnalogSample:(const void*)daqCount byteOrder:(CFByteOrder)byteOrder; //output is in byteOrder format

/*!
 @method     
 @abstract   Convert ADC data to analog (floating point) data
 @discussion Converts raw ADC data, as read from the ADC device to analog (floating point) data in physical units. External device gain and physical units are defined by an instance of id<AUIExternalDevice> associated with the stream that produced inData.
 
 External gain defines the mapping from voltage to physical units. Specifically voltage/externalGain should give physical units.
 
 Output data is floating point value:
     outData[i] = (inData[i]/voltToSampleMult)/externalGain
 with outData in byteOrder endian format.
 
 @param inData Raw data (array of unsigned integer samples) read from DAC device. Data must be native-endian format. type by ANALOG_USER_SAMPLE_TYPE.
 @param sampleBytes Number of bytes per sample of raw data (e.g. 2 for unsigned short samples)
 @param externalGain External gain, as returned from id<AUIExternalDevice>.inputGain.
 @param voltToSampleMult Scalar that converts DAC samples voltage.
 @param byteOrder Output byte order (big/little endian)
 
 @result NSData wrapping ANALOG_USER_SAMPLE_TYPE array of floating point samples in byteOrder endian format with physical units.
 
 @updated 12-10-2007
*/

+ (NSData*)daqToAnalogData:(NSData*)inData 
               sampleBytes:(unsigned)sampleBytes 
        externalDeviceGain:(double)externalGain 
          voltToSampleMult:(double)voltToSampleMult 
                 byteOrder:(CFByteOrder)byteOrder;


/*!
 @method     
 @abstract   Convert floating point analog data to DAC samples for output.
 @discussion Converts floating point analog data (in physical units) to DAC samples. External device gain and physical units are defined by an instance of id<AUIExternalDevice> associated with the stream that produced inData.
 
 External gain defines the mapping from voltage to physical units. Specifically physcialUnits/externalGain should give voltage.
 
 Output data is integer:
 outData[i] = (inData[i]/externalGain)*voltToSampleMult;
 with outData in native endian format.
 
 @param analogData Native endian floating point data (in physical units). type by ANALOG_USER_SAMPLE_TYPE.
 @param externalGain External gain, as returned from id<AUIExternalDevice>.outputGain
 @param voltToSampleMult Scalar that converts voltage to DAC samples.
 @param sampleBytes Number of bytes (unsigned integer) type of output data.
 
 @result NSData wrapping unsigned integer (of sampleBytes size) array of DAC samples.
 
 @updated 12-10-2007
 */

+ (NSData*)analogDataToDAQ:(NSData*)analogData
        externalDeviceGain:(double)externalGain
          voltToSampleMult:(double)voltToSampleMult
               sampleBytes:(unsigned)sampleBytes ;


/*
 these will convert data between Volts and daq samples, if necessary:
 1. divide by external gain    
 2. multiply by DAQ controllers sampleToVoltMult or voltToSampleMult
 3. Convert to sample (sampleBytes)
 */
-(NSData*)pullBytes:(unsigned)length byteOrder:(CFByteOrder)byteOrder; //must be input channel. data contains ANALOG_USER_SAMPLE_TYPE wide floats for analog, [controller sampleBytes] wide integers for digital. saves last sample as last sample, if this call will empty the buffer. data is passed through externalDevice's gain before being returned. Units are in externalDevice's unitsForInput.
-(void)pushBytes:(NSData*)d; //Data must be ANALOG_USER_SAMPLE_TYPE wide floats for analog, [controller sampleBytes] wide integers for digital. must be output channel. data is passed through external devices outputGain. Units passed in should be in externalDevice's unitsForOutput.

-(NSNumber*)pullLastSample; //this is probably what you want to use to get lastSample: flushes data if it's not empty and returns self.lastSample
-(void)flushDataSaveLastSample:(BOOL)save; //saves last sample as self.lastSample
-(void)flushDataBytes:(unsigned long)nBytes; //flushes up to nBytes. saves last sample if nBytes>= current buffer length

//- (BOOL)isDAQStreamInfoInitialized;
//- (void)setDAQStreamInfoInitialized:(BOOL)flag;
@end