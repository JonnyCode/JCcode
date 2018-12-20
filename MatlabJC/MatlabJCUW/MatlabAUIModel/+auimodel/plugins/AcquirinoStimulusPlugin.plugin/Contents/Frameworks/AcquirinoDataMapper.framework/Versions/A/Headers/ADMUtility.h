//
//  ADMUtility.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/20/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>
#import <Acquirino/ITCAcq.h>

@interface ADMUtility : NSObject {

}

// some parameter structures have 12 bytes of uninitialized (spare) space at
// the end.  this obviously makes it impossible to exactly check the equality
// of params, but the first chunk of bytes being equal is a pretty good
// indicator for testing purposes.  this method compares the first chunk of
// bytes, which in some cases are the only significant ones.
+ (BOOL) stimParams: (NSData *) dataA equalWithinTolerance: (NSData *) dataB;

+ (NSPersistentStoreCoordinator *) newInMemoryAUICoordinator;
+ (NSManagedObjectContext *) newInMemoryAUIContext;

@end

