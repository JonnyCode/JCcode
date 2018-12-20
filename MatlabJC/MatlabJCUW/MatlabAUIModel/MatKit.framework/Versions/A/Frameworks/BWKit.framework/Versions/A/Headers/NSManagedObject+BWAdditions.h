//
//  NSManagedObject+BWAdditions.h
//  BWKit
//
//  Created by Barry Wark on 9/23/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSManagedObject (BWAdditions)
/*!
    @method     
    @abstract   Refresh this object in it's managed object context.
    @discussion Wrapper for [[self managedObjectContext] refreshObject:self mergeChanges:merge].
*/

- (void)refreshMergeChanges:(BOOL)merge;
@end
