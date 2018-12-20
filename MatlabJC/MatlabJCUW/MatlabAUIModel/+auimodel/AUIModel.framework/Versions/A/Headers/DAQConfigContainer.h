//
//  DAQConfigContainer.h
//  AcqUI
//
//  Created by Barry Wark on 9/5/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Epoch;

@interface DAQConfigContainer :  NSManagedObject  
{
}

- (NSData *)daqConfigData;
- (void)setDaqConfigData:(NSData *)value;

// Access to-many relationship via -[NSObject mutableSetValueForKey:]
- (void)addEpochsObject:(Epoch *)value;
- (void)removeEpochsObject:(Epoch *)value;

@end
