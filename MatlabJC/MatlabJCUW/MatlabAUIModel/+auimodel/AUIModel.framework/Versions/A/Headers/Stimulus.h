//
//  Stimulus.h
//  AcqUI
//
//  Created by Barry Wark on 12/11/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import "IOBase.h"

@class Epoch;

@interface Stimulus :  IOBase  
{
}

//- (UNKNOWN_TYPE)parameters;
//- (void)setParameters:(UNKNOWN_TYPE)value;

- (NSNumber *)sampleRate;
- (void)setSampleRate:(NSNumber *)value;

- (NSNumber *)duration;
- (void)setDuration:(NSNumber *)value;

- (NSData *)parametersData;
- (void)setParametersData:(NSData *)value;

- (NSNumber *)version;
- (void)setVersion:(NSNumber *)value;

- (NSString *)stimulusID;
- (void)setStimulusID:(NSString *)value;

- (NSString *)stimDescription;
- (void)setStimDescription:(NSString *)value;

- (NSData *)data;
- (void)setData:(NSData *)value;

- (Stimulus *)parentStimulus;
- (void)setParentStimulus:(Stimulus *)value;

- (Epoch *)epoch;
- (void)setEpoch:(Epoch *)value;

@property (retain) NSSet* subStimuli;

@end

// coalesce these into one @interface Stimulus (CoreDataGeneratedAccessors) section
@interface Stimulus (CoreDataGeneratedAccessors)
- (void)addSubStimuliObject:(Stimulus *)value;
- (void)removeSubStimuliObject:(Stimulus *)value;
- (void)addSubStimuli:(NSSet *)value;
- (void)removeSubStimuli:(NSSet *)value;

@end

