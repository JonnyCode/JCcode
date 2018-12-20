//
//  NSArray-NSSetAdditions.h
//  VizUI
//
//  Created by Barry Wark on 8/30/07.
//  Copyright 2007 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSArray (NSSetAdditions)
+ (NSArray*)arrayWithSet:(NSSet*)set;
@end
