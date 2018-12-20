//
//  BWFileSystemResource.h
//  BWKit
//
//  Created by Barry Wark on 10/20/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>


@interface BWFileSystemResource :  NSManagedObject  
{
}

@property (retain,nonatomic) id url;
@property (retain) id alias;
@property (readonly,nonatomic) NSString *path;
@property (retain) NSURL *relativeURL;

/**
 Root for relative Alias resolution.
 
 Returns nil by default (ie using an Alias that does not start its resolution
 search relative to a particular root). Subclasses may override this method to
 Return a relative root.
 
 @return Path to search root.
 */
- (NSURL*)relativeRoot;


/**
 Uses provided relativeRoot. For migration.
 */
+ (void)setURL:(NSURL*)url relativeToRootURL:(NSURL*)relativeRoot forFileSystemResource:(NSManagedObject*)resource;
+ (NSURL*)URLForFileSystemResource:(NSManagedObject*)resource relativeToRootURL:(NSURL*)relativeRoot;
@end


