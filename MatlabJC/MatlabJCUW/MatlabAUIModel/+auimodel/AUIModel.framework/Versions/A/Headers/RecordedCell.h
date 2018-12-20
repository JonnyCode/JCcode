//
//  RecordedCell.h
//  AcqUI
//
//  Created by Barry Wark on 12/11/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Epoch;
@class Experiment;

@interface RecordedCell :  NSManagedObject  
{
}

- (NSString *)label;
- (void)setLabel:(NSString *)value;

- (NSDate *)startDate;
- (void)setStartDate:(NSDate *)value;

- (NSString *)comment;
- (void)setComment:(NSString *)value;

// Access to-many relationship via -[NSObject mutableSetValueForKey:]
- (void)addEpochsObject:(Epoch *)value;
- (void)removeEpochsObject:(Epoch *)value;

- (Experiment *)experiment;
- (void)setExperiment:(Experiment *)value;

@end
