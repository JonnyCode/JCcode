//
//  NSArray+SortAdditions.h
//  VizUI
//
//  Created by Barry Wark on 10/11/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSArray (SortAdditions)
- (NSArray*)sortedArrayUsingKey:(NSString*)key ascending:(BOOL)ascending;
@end
