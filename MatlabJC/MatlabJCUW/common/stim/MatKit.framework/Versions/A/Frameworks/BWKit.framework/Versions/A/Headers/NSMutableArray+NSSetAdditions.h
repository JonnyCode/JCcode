//
//  NSMutableArray+NSSetAdditions.h
//  VizUI
//
//  Created by Barry Wark on 10/1/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSMutableArray (NSSetAdditions)
+ (NSMutableArray*)arrayWithSet:(NSSet*)set;
@end
