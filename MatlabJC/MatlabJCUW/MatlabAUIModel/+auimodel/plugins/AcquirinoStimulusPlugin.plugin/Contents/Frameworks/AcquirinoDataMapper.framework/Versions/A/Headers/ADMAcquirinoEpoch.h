//
//  ADMAcquirinoEpoch.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/20/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <BWKit/BWKit.h>

// simplified, decoupled constants for amplifier mode.  exhaustive version is in
// AcqUI/AUIAxoPatchDevice.h.  values correspond to those in Acquirino data file
typedef enum {
  CURRENT = 2,
  VOLTAGE = 6
} ADMAmpMode;

@interface ADMAcquirinoChannel : NSObject {
@private
  char segmentNumber;
  NSNumber* channelNumber;
}

@property (readwrite) char segmentNumber;
@property (readwrite,assign) NSNumber* channelNumber;

@end

// data objects for graph returned by the enumerator
@interface ADMAcquirinoInputChannel : ADMAcquirinoChannel {
@private
  ADMAmpMode ampMode;
  NSNumber* ampGain;
  NSData* sampledAdcCounts;
}

@property (readwrite,assign) ADMAmpMode ampMode;
@property (readwrite,assign) NSNumber* ampGain;
@property (readwrite,assign) NSData* sampledAdcCounts;

@end

@interface ADMAcquirinoOutputChannel : ADMAcquirinoChannel {
@private
  NSNumber* stimulusTypeCode;
  NSData* stimulusParameters;
};

@property (readwrite,assign) NSNumber* stimulusTypeCode;
@property (readwrite,assign) NSData* stimulusParameters;
@end

@interface ADMAcquirinoGClampOutputChannel : ADMAcquirinoOutputChannel {
  BWNumericData* excitatoryWaveform;
  BWNumericData* inhibitoryWaveform;
};

@property (readwrite,assign) BWNumericData* excitatoryWaveform;
@property (readwrite,assign) BWNumericData* inhibitoryWaveform;
@end

@interface ADMAcquirinoEpoch : NSObject {
@private
  NSNumber* number;
  NSNumber* time;
  NSString* comment;
  NSNumber* samplingInterval;
  NSNumber* dataPointCount;
  NSDictionary* outputChannels; 
  NSDictionary* inputChannels; 
}

@property (readwrite,assign) NSNumber* number;
@property (readwrite,assign) NSNumber* time;
@property (readwrite,assign) NSString* comment;
@property (readwrite,assign) NSNumber* samplingInterval;
@property (readwrite,assign) NSNumber* dataPointCount;
@property (readwrite,assign) NSDictionary* outputChannels;
@property (readwrite,assign) NSDictionary* inputChannels;

@end

