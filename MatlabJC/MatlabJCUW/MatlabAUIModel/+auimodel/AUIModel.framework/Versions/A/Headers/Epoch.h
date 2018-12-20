//
//  Epoch.h
//  AcqUI
//
//  Created by Barry Wark on 11/16/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Response;
@class Stimulus;
@class Resource;
@class DAQConfigContainer;
@class RecordedCell;
@class KeywordTag;

@interface Epoch :  NSManagedObject  
{
}

@property (retain) NSNumber * userDouble3;
@property (retain) NSNumber * duration;
@property (retain) NSString * protocolID;
@property (retain) NSNumber * userDouble5;
@property (retain) NSDate * startDate;
@property (retain) NSNumber * userDouble1;
@property (retain) NSNumber * includeInAnalysis;
@property (retain) NSNumber * saveResponse;
@property (retain) NSNumber * userDouble2;
@property (retain) NSString * userString;
@property (retain) NSNumber * bathTemperature;
@property (retain) NSDictionary* protocolSettings;
@property (retain) NSString * comment;
@property (retain) NSNumber * userDouble4;
@property (retain) NSSet* responses;
@property (retain) NSSet* stimuli;
@property (retain) NSSet* protocolSettingsKVPairs;
@property (retain) NSSet* resources;
@property (retain) DAQConfigContainer * daqConfig;
@property (retain) RecordedCell * cell;
@property (retain) NSSet* keywords;

@end

@interface Epoch (CoreDataGeneratedAccessors)
- (void)addResponsesObject:(Response *)value;
- (void)removeResponsesObject:(Response *)value;
- (void)addResponses:(NSSet *)value;
- (void)removeResponses:(NSSet *)value;

- (void)addStimuliObject:(Stimulus *)value;
- (void)removeStimuliObject:(Stimulus *)value;
- (void)addStimuli:(NSSet *)value;
- (void)removeStimuli:(NSSet *)value;

- (void)addProtocolSettingsKVPairsObject:(NSManagedObject *)value;
- (void)removeProtocolSettingsKVPairsObject:(NSManagedObject *)value;
- (void)addProtocolSettingsKVPairs:(NSSet *)value;
- (void)removeProtocolSettingsKVPairs:(NSSet *)value;

- (void)addResourcesObject:(Resource *)value;
- (void)removeResourcesObject:(Resource *)value;
- (void)addResources:(NSSet *)value;
- (void)removeResources:(NSSet *)value;

- (void)addKeywordsObject:(KeywordTag *)value;
- (void)removeKeywordsObject:(KeywordTag *)value;
- (void)addKeywords:(NSSet *)value;
- (void)removeKeywords:(NSSet *)value;
@end

