//
//  Experiment.h
//  AcqUI
//
//  Created by Barry Wark on 9/6/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class RecordedCell;
@class Note;
@class AUIFileSystemResource;
@class KeywordTag;

@interface Experiment : NSManagedObject
{
}

@property (retain) NSString * daqID;
@property (retain) NSString * otherNotes;
@property (retain) NSString * purpose;
@property (retain) NSData * rigSettingsData;
@property (retain) NSDate * startDate;
@property (retain) NSSet* cells;
@property (retain) NSSet* keywords;
@property (retain) NSSet* notes;
@property (retain) NSSet* resources;
@property (retain) AUIFileSystemResource * responseDataFile;

/**
 Create a response data file at given URL.
 
 @param url URL for response data file. Must be a fileURL. @see responseDataURLForPersistentStoreURL:.
 */
+ (void)createResponseDataFileAtURL:(NSURL*)url;

/**
 Calculate the response data file URL for a given persistent store's Experiment.
 
 Calculates an unused file name in the desired location. There is therefore a race condition if pStoreURL is in a public directory.
 
 @param pStoreURL URL of an Experiment's persistent store.
 @return URL for the resposne data file.
 */
+ (NSURL*)responseDataURLForPersistentStoreURL:(NSURL*)pStoreURL;

/**
 Use a given URL for the response data file.
 Creates a file system resource referencing the given URL and sets that resource as self.responseDataFile.
 
 @throw NSAssertionError if there's already a self.responseDataFile.
 @param url Must be fileURL to an existing file.
 */
- (void)useResponseDataFileAtURL:(NSURL*)url;
@end

// coalesce these into one @interface Experiment (CoreDataGeneratedAccessors) section
@interface Experiment (CoreDataGeneratedAccessors)
- (void)addCellsObject:(RecordedCell *)value;
- (void)removeCellsObject:(RecordedCell *)value;
- (void)addCells:(NSSet *)value;
- (void)removeCells:(NSSet *)value;

- (void)addKeywordsObject:(KeywordTag *)value;
- (void)removeKeywordsObject:(KeywordTag *)value;
- (void)addKeywords:(NSSet *)value;
- (void)removeKeywords:(NSSet *)value;

- (void)addNotesObject:(Note *)value;
- (void)removeNotesObject:(Note *)value;
- (void)addNotes:(NSSet *)value;
- (void)removeNotes:(NSSet *)value;

- (void)addResourcesObject:(NSManagedObject *)value;
- (void)removeResourcesObject:(NSManagedObject *)value;
- (void)addResources:(NSSet *)value;
- (void)removeResources:(NSSet *)value;

@end
