//
//  ADMStimulusParameterMapper.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/17/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

// constructed from constants #defined in ITCAcq.h, VERSIONNUMBER = 8
#define ADM_MAX_STIMULUS_CODE 13

extern NSString* const ADMAcquirinoStimulusTypeKey;
extern NSString* const ADMAcquirinoStimulusTypeNameKey;

extern NSString * const ADMAcquirinoGClampStimulusConductanceWaveFormKey;
extern NSString * const ADMAcquirinoGClampStimulusReversalPotentialKey;
extern NSString * const ADMAcquirinoGClampStimulusConductanceTypeKey;

extern int const stimulusCode2ParametersSize[ADM_MAX_STIMULUS_CODE + 1];
extern NSString* const stimulusCode2String[ADM_MAX_STIMULUS_CODE + 1];

@interface ADMStimulusParameterMapper : NSObject {
}

+ (NSMutableDictionary *) dictionaryForStimulusType: (NSNumber *) type
  andParameters: (NSData *) parameters;

+ (NSData *) parametersForDictionary: (NSDictionary *) pairs
  andStimulusType: (NSNumber *) type;

@end

