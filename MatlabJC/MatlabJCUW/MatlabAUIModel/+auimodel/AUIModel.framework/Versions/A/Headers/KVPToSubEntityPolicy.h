//
//  KVPToSubEntityPolicy.h
//  AcqUI
//
//  Created by Barry Wark on 12/17/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface KVPToSubEntityPolicy : NSEntityMigrationPolicy {

}

- (BOOL)createDestinationInstancesForSourceInstance:(NSManagedObject *)sInstance 
                                      entityMapping:(NSEntityMapping *)mapping 
                                            manager:(NSMigrationManager *)manager 
                                              error:(NSError **)error;

- (BOOL)performCustomValidationForEntityMapping:(NSEntityMapping *)mapping 
                                        manager:(NSMigrationManager *)manager 
                                          error:(NSError **)error;

@end
