//
//  Project.h
//  Ovation
//
//  Created by Barry Wark on 6/9/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class DataStore;

@interface Project :  NSManagedObject  
{
}

@property (retain) id uiPreferences;
@property (retain) id projectRootURL;
@property (retain) NSString * projectDescription;
@property (retain) NSPredicate * searchPredicate;
@property (retain) NSSet* dataStores;
@property (retain) NSSet* savedPredicates;
@property (retain) id userInfo;
@end

@interface Project (CoreDataGeneratedAccessors)
- (void)addDataStoresObject:(DataStore *)value;
- (void)removeDataStoresObject:(DataStore *)value;
- (void)addDataStores:(NSSet *)value;
- (void)removeDataStores:(NSSet *)value;

- (void)addSavedPredicatesObject:(NSManagedObject *)value;
- (void)removeSavedPredicatesObject:(NSManagedObject *)value;
- (void)addSavedPredicates:(NSSet *)value;
- (void)removeSavedPredicates:(NSSet *)value;

@end

