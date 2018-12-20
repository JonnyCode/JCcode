//
//  DataStore.h
//  Ovation
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Project;

@interface DataStore :  NSManagedObject  
{
}

@property (retain,readonly) NSSet* persistentStoreURLs;
@property (retain) Project * project;


- (BOOL)treeIsLeaf;
@end


