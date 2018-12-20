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
 
 Entity userInfo may also define a relativeRootURIKeypath entry giving the key
 path to use as the relative root. If defined, this value of this keypath is used
 instead of nil.
 
 @return Path to search root.
 */
- (NSURL*)relativeRoot;


/**
 Uses provided relativeRoot. For migration.
 */
+ (void)setURL:(NSURL*)url relativeToRootURL:(NSURL*)relativeRoot forFileSystemResource:(BWFileSystemResource*)resource;
+ (NSURL*)URLForFileSystemResource:(BWFileSystemResource*)resource relativeToRootURL:(NSURL*)relativeRoot;
+ (NSURL*)URLForFileSystemResource:(BWFileSystemResource*)resource relativeToRootURL:(NSURL*)relativeRoot openPanel:(NSOpenPanel*)op;
@end


