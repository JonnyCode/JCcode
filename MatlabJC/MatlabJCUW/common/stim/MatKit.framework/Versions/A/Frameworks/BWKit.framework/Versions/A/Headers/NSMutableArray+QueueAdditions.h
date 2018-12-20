//
//  NSMutableArray+QueueAdditions.h
//  BWKit
//
//  Created by Barry Wark on 1/18/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSMutableArray (QueueAdditions)
/*!
    @method     
    @abstract   Remove and return first object
    @discussion Atomic.
*/

- (id)popFirst;


/*!
    @method     
    @abstract   Remove and return last object in array (last push'd object)
    @discussion Atomic.
*/

- (id)pop;

/*!
    @method     
    @abstract   Add an object to the end of the array (push).
    @discussion Synonym for self.addObject:. Atomic.
*/

- (void)push:(id)obj;

@end
