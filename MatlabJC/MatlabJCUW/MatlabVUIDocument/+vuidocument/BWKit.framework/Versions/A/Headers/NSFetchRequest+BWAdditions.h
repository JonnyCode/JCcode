//
//  NSFetchRequest+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 7/7/07.
//  Copyright 2007 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSFetchRequest (BWAdditions)
+ (NSFetchRequest*)fetchRequestWithPredicate:(NSPredicate*)predicate
                                 entityName:(NSString*)entityName
                       managedObjectContext:(NSManagedObjectContext*)moc;

/*!
 @method     
 @abstract   Factory method for creating NSFetchRequeste
 @discussion Depreciated. Use fetchRequestForEntityName:managedObjectContext:predicateFormat:
 @depreciated
*/

+ (NSFetchRequest*)fetchRequestForEntityName:(NSString*)entityName
                        managedObjectContext:(NSManagedObjectContext*)moc
                             predicateFormat:(NSString*)format,...;
@end
