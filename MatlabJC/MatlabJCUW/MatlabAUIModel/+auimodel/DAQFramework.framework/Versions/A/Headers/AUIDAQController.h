//
//  AUIDAQController.h
//  AcqUI
//
//  Created by Barry Wark on 3/31/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "DAQFramework/AUIDAQStream.h"
#import "DAQFramework/AUIIOController.h"

#define DAQ_DEFAULT_SAMPLE_RATE 10000.0
#define DAQ_CONTROLLER_EXCEPTION_NAME @"DAQ_CONTROLLER_EXCEPTION_NAME"


@interface AUIDAQController : AUIIOController { //abstract base class for DAQControllers 
	//double sampleRate;
    NSMutableSet *channelsInUse; //channels that are will be serviced by readFromFIFO
    
	BOOL hardwareInitialized;
	BOOL isRunning;
}

/*!
    @method     
    @abstract   Initialize channels and (optionally) initHardware.
    @discussion Initializes channels set and channelsInUse set. Subclasses should add channels to the channels set as appropriate.
*/

- (id)initWithNumChannels:(unsigned)nChannels initHardware:(BOOL)initHardware;

//append data to be written to channel. throws AUIIOException if channel is not available. Data must be host-native endian format. See AUIDAQStreamInfo for other requirements.
-(void)writeToFIFO:(unsigned short)channel ofType:(StreamType)type data:(NSData*)data;

//read nsamples off channel. runs IO until nsamples are available on channel (appends data onto output channel until read).
// rChannels is list of all channels to read from. Data is in Little endian format
-(NSData*)readFromFIFO:(unsigned short)channel ofType:(StreamType)type numSamples:(unsigned long)nsamples readStreams:(NSSet*)rStreams;

//all are in DAQ counts // !!!:barry:20070709 todo: these need to be updated to use correct units, with DAQStreamInfo providing conversion
-(void)asyncWriteToDAC:(unsigned short)channel value:(short)v; //deprecated 0.4
-(void)asyncWriteToDigital:(unsigned short)channel value:(short)v; //deprecated 0.4
-(void)asyncWriteToAux:(unsigned short)channel value:(short)v; //deprecated 0.4
-(NSNumber*)asyncReadFromADC:(unsigned short)channel; //deprecated 0.4
-(NSNumber*)asyncReadFromDigital:(unsigned short)channel; //deprecated 0.4
-(NSNumber*)asyncReadFromAux:(unsigned short)channel; //deprecated 0.4

- (void)asyncWriteSample:(NSNumber*)sample
                toStream:(AUIDAQStream*)s;
- (NSNumber*)asyncReadFromStream:(AUIDAQStream*)s;

-(ANALOG_USER_SAMPLE_TYPE)asyncReadSampleFromAnalogStream:(AUIDAQStream*)stream; //returns value in device physical units

// read from aux (digital) channel. returns a dictionary with provided keys and aux values as return (bool in NSNumber)
-(NSDictionary*)asyncReadFromAux:(unsigned short)channel keys:(NSArray*)keys;
//similarly, write to aux (digital) channel. values of d must be NSNumber (should be bool value)
-(void)asyncWriteToAux:(unsigned short)channel keys:(NSDictionary*)d;
-(short)bitForAuxKey:(id)k; //returns -1 if not found

- (double)minimumSamplingRateHz;
- (double)maximumSamplingRateHz;

-(unsigned long)availableFIFOSamples:(unsigned short)channel ofType:(StreamType)type; //number remaining in write or number available in read
-(unsigned long)availableStreamSamples:(unsigned short)channel ofType:(StreamType)type; //number in StreamInfo NSData buffer

//- (double)sampleRate;
//- (void)setSampleRate:(double)newSampleRate;

- (NSMutableSet *)channelsInUse;
- (void)setChannelsInUse:(NSMutableSet *)newChannelsInUse;
- (void)addToChannelsInUse:(id)channelsInUseObject;
- (void)removeFromChannelsInUse:(id)channelsInUseObject;

+(id)pluginType;
@end
