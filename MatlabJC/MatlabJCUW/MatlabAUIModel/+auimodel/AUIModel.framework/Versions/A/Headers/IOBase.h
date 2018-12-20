//
//  IOBase.h
//  AcqUI
//
//  Created by Barry Wark on 12/17/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Resource;

@interface IOBase :  NSManagedObject  
{
}


@property (retain) NSNumber * externalDeviceMode;
@property (retain) NSNumber * type;
@property (retain) NSNumber * externalDeviceGain;
@property (retain) NSNumber * channelID;
@property (retain) NSNumber * sampleBytes;
@property (retain) NSSet* resources;
@property (retain) NSSet* keywords;

@end

@interface IOBase (CoreDataGeneratedAccessors)
- (void)addResourcesObject:(Resource *)value;
- (void)removeResourcesObject:(Resource *)value;
- (void)addResources:(NSSet *)value;
- (void)removeResources:(NSSet *)value;

- (void)addKeywordsObject:(NSManagedObject *)value;
- (void)removeKeywordsObject:(NSManagedObject *)value;
- (void)addKeywords:(NSSet *)value;
- (void)removeKeywords:(NSSet *)value;

@end

