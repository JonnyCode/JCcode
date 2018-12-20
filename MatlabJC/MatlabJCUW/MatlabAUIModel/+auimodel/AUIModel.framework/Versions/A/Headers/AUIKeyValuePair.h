//
//  AUIKeyValuePair.h
//  AcqUI
//
//  Created by Barry Wark on 11/7/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/KeyValuePair.h>

@class Stimulus;
@class Epoch;

@interface AUIKeyValuePair : KeyValuePair
{
}

@property (retain) NSSet* epochs;
@property (retain) NSSet* stimuli;
@property (retain) NSManagedObject * streamValue;

@end

// coalesce these into one @interface AUIKeyValuePair (CoreDataGeneratedAccessors) section
@interface AUIKeyValuePair (CoreDataGeneratedAccessors)
- (void)addEpochsObject:(Epoch *)value;
- (void)removeEpochsObject:(Epoch *)value;
- (void)addEpochs:(NSSet *)value;
- (void)removeEpochs:(NSSet *)value;

- (void)addStimuliObject:(Stimulus *)value;
- (void)removeStimuliObject:(Stimulus *)value;
- (void)addStimuli:(NSSet *)value;
- (void)removeStimuli:(NSSet *)value;

@end
