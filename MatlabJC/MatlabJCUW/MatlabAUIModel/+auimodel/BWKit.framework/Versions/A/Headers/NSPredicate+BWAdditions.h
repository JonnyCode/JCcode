//
//  NSPredicate+BWAdditions.h
//  BWKit
//
//  Created by Barry Wark on 8/2/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSPredicate (BWAdditions)
- (BOOL)isCompoundPredicate;
- (BOOL)isComparisonPredicate;
@end
