//
//  NSArray+IndexingAdditions.h
//  BWKit
//
//  Created by Barry Wark on 10/31/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSArray (IndexingAdditions)

/*!
    @method     
    @abstract   Returns object at index 0 in array (if it exists)
    @discussion Shorthand for objectAtIndex:0. Returns nil if there are no objects in the array.
*/

- (id)firstObject;
@end
