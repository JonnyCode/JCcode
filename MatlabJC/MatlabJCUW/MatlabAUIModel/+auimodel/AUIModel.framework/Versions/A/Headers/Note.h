//
//  Note.h
//  AcqUI
//
//  Created by Barry Wark on 2/14/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Experiment;

@interface Note :  NSManagedObject  
{
}

- (NSDate *)date;
- (void)setDate:(NSDate *)value;

- (NSString *)text;
- (void)setText:(NSString *)value;

- (Experiment *)experiment;
- (void)setExperiment:(Experiment *)value;

@end
