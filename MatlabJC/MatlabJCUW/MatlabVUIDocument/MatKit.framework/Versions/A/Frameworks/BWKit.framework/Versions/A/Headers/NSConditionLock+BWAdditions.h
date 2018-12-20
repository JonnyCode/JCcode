//
//  NSConditionLock-SynchronizeAdditions.h
//  AcqUI
//
//  Created by Barry Wark on 1/5/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSConditionLock (BWAdditions)
- (void)synchronizeWithCondition:(NSInteger)condition;
- (BOOL)synchronizeWithCondition:(NSInteger)condition beforeDate:(NSDate*)d;
- (void)synchronizeSetCondition:(NSInteger)condition;
- (void)synchronizeIncrementCondition;
- (void)synchronizeDecrementCondition;
@end
