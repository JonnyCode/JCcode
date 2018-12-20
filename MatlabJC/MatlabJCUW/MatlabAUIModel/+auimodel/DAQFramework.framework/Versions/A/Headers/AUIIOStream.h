//
//  AUIIOStreamInfo.h
//  AcqUI
//
//  Created by CJ Bell on 7/10/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIExternalDevice.h>

@class AUIIOController;

extern NSString *AUIIOStreamChangedNotification;

// Stream types
typedef enum {
	UNCONFIGURED=0,
	ANALOG_IN,
	ANALOG_OUT,
	DIGITAL_IN,
	DIGITAL_OUT,
    AUX_IN,
    AUX_OUT,
    VIDEO_IN,
    VIDEO_OUT
} StreamType;

typedef enum {
    UNKNOWN_TYPE=0,
    INPUT_TYPE,
    OUTPUT_TYPE
} StreamMetaType;

// Status values for OUT streams
typedef enum {
	IDLE=0,
	WRITING, // last sample has not been sent to FIFO
	LAST_SAMPLE_WRITTEN, // last sample has been sent to FIFO
	TERMINATED // background sample sent to FIFO, DAC will not update after this
} StreamStatus;

/*!
 @class
 @abstract    Base class for I/O streams in the DAQFrameworkd
 @discussion  Subclasses must call [self postStreamChangedNotification] on any changes that affect I/O from this stream
 (e.g. sampling rate).
*/

@interface AUIIOStream : NSObject <NSCopying> {
    StreamType	type;
	unsigned long	channelNumber; 		// type and number are used to map the StreamInfo to a hardware channel
	NSString *userDescription;
    BOOL enabled;
    
    AUIIOController *controller; //not retained
    
    id<AUIExternalDevice> externalDevice;
    id<AUIExternalDevice> externalDeviceOwner; //non-retained. if non-nil, this stream is acting as a secondary input for an external device (e.g. amp mode) and should not be considered as an available stream for I/O. AUIDocument should add all owned input streams to the list of read streams when acquiring. // !!!:barry:20070709 check that this is followed for AUIDocument and external devices. DAQConfig dialog should also disable device popup for these streams
    
    NSMutableDictionary *userInfo;
}

// !!!:barry:20071224 readonly so that hash cannot change after initialization.
@property (assign,readonly) StreamType type;
@property (assign,readonly) unsigned long channelNumber;
@property (assign,readonly) AUIIOController *controller; //not retained

@property (copy,readwrite) NSString *userDescription;

/*!
 @abstract userInfo dictionary
 @discussion NOT stored in stream properties.
 */
@property (retain,readwrite) NSMutableDictionary *userInfo;

/*!
 @method     
 @abstract   Creates an unbound (does not belong to a controller) AUIIOStreamInfo
 @discussion Intended to be used when the controller is not available for de-serializing a stream info.
 @param type StreamType
 @param channelNumber Channel number
 @param properties Dictionary of stream properties. If non-nil, unbound stream will have its properties set from this dictionary.
 */
+ (AUIIOStream*)unboundStreamInfoWithType:(StreamType)theType
                            channelNumber:(unsigned long)theChannelNumber
                               properties:(NSDictionary*)properties;

/*!
 @method     
 @abstract   Convenience method for unbound stream with nil properties
 @discussion Calls unboundStreamInfoWithType:channelNumber:properties: with nil properties.
*/

+ (AUIIOStream*)unboundStreamInfoWithType:(StreamType)theType
                            channelNumber:(unsigned long)theChannelNumber;

+ (BOOL)isInputStreamType:(StreamType)type;
+ (BOOL)isOutputStreamType:(StreamType)type;
+ (StreamMetaType)streamMetaType:(StreamType)type;

/*!
    @method     
    @abstract   Designated initializer for AUIIOStream
    @discussion channelNumber and type are immutable (since they're used for hasing and equality testing)
*/

- (id)initWithIOController:(AUIIOController*)c
             channelNumber:(unsigned long)theChannelNumber
                      type:(StreamType)theType; //init with default values. must be set by controller

/*
 Should return an array of KVO keys corresponding to the properties of this stream that the controller should archive.
 Controller will remove channelNumber, type, and controller automatically, so they don't need to be included.
 Subclasses may override streamPropertyKeys to add additional keys, but should append these added
 keys to the superclass' streamPropertyKeys.
 */
- (NSArray*)streamPropertyKeys;

/**
Returns YES if this stream is equal to other stream (type, channelNumber and controllerID are equal)
 */
- (BOOL)isEqualToStream:(AUIIOStream*)otherStream;
- (NSComparisonResult)compare:(AUIIOStream *)otherStream;
- (unsigned)hash;
- (BOOL)isEqual:(id)anObject;

+ (NSArray*)streamTypeNames;
- (NSString*)typeName;
+ (NSString*)typeNameForType:(StreamType)theType;
- (BOOL)isInputChannel;
- (BOOL)isOutputChannel;
- (BOOL)isAnalogChannel;
- (BOOL)isVideoChannel;
- (BOOL)isConfigured;

- (double)samplingRate;
- (void)setSamplingRate:(double)newSamplingRate;
- (BOOL)validateSamplingRate:(id *)ioValue error:(NSError **)outError;

- (NSString *)userDescription;
- (void)setUserDescription:(NSString *)anUserDescription;

- (BOOL)enabled;
- (void)setEnabled:(BOOL)flag;

- (id<AUIExternalDevice>)externalDeviceOwner;
- (void)setExternalDeviceOwner:(id<AUIExternalDevice>)anExternalDeviceOwner;

- (id<AUIExternalDevice>)externalDevice;
- (void)setExternalDevice:(id<AUIExternalDevice>)anExternalDevice;

- (NSString*)externalDeviceUnits;

- (void)flushData;
- (void)flushDataSaveLastSample:(BOOL)save;

- (NSData*)backgroundDataForSamples:(unsigned long)nSamples;

/*!
 An archiveable representation of this stream's properties
 */
- (id)streamProperties;

/*!
 An archiveable identifier that uniquely ID's this stream (encapsulates channelNumber, type, controllerID)
 */
- (id)identifier;

/*!
    @method     
    @abstract   Alias for self.identifier()
    @discussion An archiveable identifier that uniquely IDs this stream (controllerID, channelNumber, type)
*/

- (id)streamIdentifier;

/*!
    @method     
    @abstract   Convenience method to return self.controller.controllerID
*/

- (id)controllerID;

/*!
    @method     
    @abstract   Test if stream is bound to a controller
    @discussion isUnbound := (self.controller == nil)
*/

- (BOOL)isUnbound;
@end
