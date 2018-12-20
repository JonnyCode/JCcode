//
//  AUIModelMigrator.h
//  AcqUI
//
//  Created by Barry Wark on 12/28/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface AUIModelMigrator : NSObject {

}

+ (AUIModelMigrator*)sharedMigrator;

- (BOOL)migrateStoreAtURL:(NSURL*)url
                    error:(NSError**)err;
- (BOOL)migrateStoreAtURL:(NSURL*)url
                    error:(NSError**)err
             showProgress:(BOOL)showProgress;
@end
