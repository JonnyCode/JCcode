//
//  KeywordTag.h
//  AcqUI
//
//  Created by Barry Wark on 1/29/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>


@interface KeywordTag :  NSManagedObject  
{
}

@property (retain) NSString * tag;

/**
 Creates (if needed) a keyword tag with given tag in owner's managed object context.
 
 If owner.objectID.persistentStore is non-nil, returned KeywordTag is unique in owner's persistent store
 and is assigned to the same persistent store.
 
 @param tag Key word tag.
 @param owner Object that will be associated with this tag.
 @param error OUT error value.
 
 @return Unique (in owner's persistent store, if non-nil) KeywordTag.
 */

+ (KeywordTag*)keywordTagWithTag:(NSString*)tag
                        forOwner:(NSManagedObject*)owner
                           error:(NSError**)error;
/**
 Creates (if needed) a unique keyword tag with given tag in owner's managed object context.
 */
+ (KeywordTag*)keywordTagWithTag:(NSString*)tag
           inManagedObjectContext:(NSManagedObjectContext*)moc
                            error:(NSError**)error;

@end


