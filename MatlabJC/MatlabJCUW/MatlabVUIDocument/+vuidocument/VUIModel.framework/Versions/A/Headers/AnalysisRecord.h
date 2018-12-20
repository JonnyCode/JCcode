//
//  AnalysisRecord.h
//  Ovation
//
//  Created by Barry Wark on 10/20/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>

@class Project;

/**
 Implementation for Ovation's VUIDocument.xcdatamodel AnalysisRecord.
 An analysis record encasulates information (data, code, parameters) of an analysis.
 */

@interface AnalysisRecord :  NSManagedObject  
{
}

@property (retain) NSString * userDescription;
@property (retain) NSString * name;
@property (retain) NSDate * timeStamp;
@property (retain) NSSet* keywords;
@property (retain) NSSet* resources;
@property (retain) NSNumber * codeRevision;
@property (retain) NSURL* svnRepositoryURI;
@property (retain) NSDictionary* parameters;
@property (readonly) NSAttributedString *attributedCodeURL;
@property (readonly) NSString *codeCheckOutCommand;
@property (retain) id resourceRootURI;


/**
 Auto-generated tag for this analylsis record. 
 
 Once initialized by the factory method, an AnalysisRecord's autoTag is an immutable property that identifies the AnalysisRecord and tags any Epochs passed into the factory method. In addition, Ovation, will assume that any relevant epochs are tagged with the AnalysisRecord's autoTag, so if they're not passed into the AnalsyisRecord factory method, you will tag them manually.
 
*/
@property (retain,readonly) NSString *autoTag;


/**
 Factory method for creating AnalysisRecord instances.
 
 Builds an AnalysisRecord instance in an NSManagedObjectContext. An 
 auto-generated keyword tag of the form <name>-<timeStammp> will be added to the
 record and resources. Epochs associated with this record will be tagged with
 this record's auto-generated tag (in the Epochs' managed object context).
 
 Although you can create an analysis record by hand, none of the automatic
 tagging, as described above, will be done for you.
 
 @param name The record name.
 @param description A (multi-line) description of the record.
 @param parameters An NSDictionary of analysis parameters. Keys must be `NSString` for display.
 @param epochs NSSet of Epoch instances (AUIModel) associated with this analysis.
 @param keywordStrings NSSet of strings for record keywords.
 @param resourceURLs NSSet of resource URLs to link to this record.
 @param resourceRoot NSURL for relative resolution of resource paths
 @param rev Code revision for this record.
 @param codeURL NSURL of code repository.
 @param moc The NSManagedObjectContext in which to insert the record.
 @param[out] err Returned error by reference or nil.
 @exception NSInternalInconsistencyException if parameter's keys are not NSString.
 @return An AnalysisRecord instance or nil if an error occurs.
 */
+ (AnalysisRecord*)insertAnalysisRecordWithName:(NSString*)name
                                    description:(NSString*)description
                                     parameters:(NSDictionary*)parameters
                                         epochs:(NSSet*)epochs
                                       keywords:(NSSet*)keywordStrings
                                      resources:(NSSet*)resourceURLs
                                resourceRootURL:(NSURL*)resourceRoot
                                   codeRevision:(NSUInteger)rev
                                        codeURL:(NSURL*)codeURL
                         inManagedObjectContext:(NSManagedObjectContext*)moc
                                          error:(NSError**)err;

/**
 * Wrapper that accepts properties in dictionary form and breaks them out into
 * parameters to the factory method.
 *
 * Each parameter to the factory method must appear as a key in the dictionary,
 * except for moc, err, and epochs.  Epochs are assumed to be tagged elsewhere,
 * and so nil is passed for this param.
 *
 * @param props Properties dict keys named as factory params.
 * @param moc VUIModel managed object context.
 * @exception NSError Any error returned by factory method.
 * @exception NSInternalInconsistencyException as thrown by factory method.
 * @return Return of factory.
 */

+ (AnalysisRecord *) insertAnalysisRecordForDictionary: (NSDictionary *) props
                              intoManagedObjectContext: (NSManagedObjectContext *) moc;

@end

@interface AnalysisRecord (CoreDataGeneratedAccessors)
- (void)addKeywordsObject:(NSManagedObject *)value;
- (void)removeKeywordsObject:(NSManagedObject *)value;
- (void)addKeywords:(NSSet *)value;
- (void)removeKeywords:(NSSet *)value;

- (void)addResourcesObject:(NSManagedObject *)value;
- (void)removeResourcesObject:(NSManagedObject *)value;
- (void)addResources:(NSSet *)value;
- (void)removeResources:(NSSet *)value;
@end
