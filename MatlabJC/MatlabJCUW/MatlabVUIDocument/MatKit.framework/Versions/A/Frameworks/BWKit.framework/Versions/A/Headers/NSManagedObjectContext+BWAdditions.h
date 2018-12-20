//
//  NSManagedObjectContext+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 5/4/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSManagedObjectContext (BWAdditions)
+ (NSManagedObjectContext*)context;
+ (NSManagedObjectContext*)contextWithPersistentStoreCoordinator:(NSPersistentStoreCoordinator*)psc;
- (NSArray*)objectsWithIDs:(NSArray*)objectIDs;

/**
 Execute a fetch request that should return at most 1 result.
 
 @return The single result or nil if no objects match request.
 @throw BWErrorException if an error is returned from the underlying executeFetchRequest:error: or if more than one object is found.
 */
- (NSManagedObject*)executeSingleResultFetch:(NSFetchRequest*)request;

- (id)objectWithURI:(NSURL*)url;
@end
