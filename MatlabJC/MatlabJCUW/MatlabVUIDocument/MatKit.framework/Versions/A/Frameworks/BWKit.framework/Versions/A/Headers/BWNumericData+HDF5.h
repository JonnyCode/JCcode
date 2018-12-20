//
//  BWNumericData+HDF5.h
//  BWKit
//
//  Created by Barry Wark on 1/30/09.
//  Copyright 2009 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWNumericData.h>

@interface BWNumericData (HDF5)
/**
 Read and return a BWNumericData instance from an HDF5 data set.
 
 @param url URL to an HDF5 file.
 @param dataSet Data set name in HDF5 file.
 @param err OUT NSError* instance if failure.
 @return BWNumericData instance of nil if error.
 */

+ (BWNumericData*)numericDataFromHDF5URL:(NSURL*)url dataSet:(NSString*)dataSet error:(NSError**)err;

/**
 Write numeric data to an HDF5 data set. Dtype information is written as a dataset sting attribute so that it can be retrieved by +numericDataFromHDF5URL:dataSet:error:.
 
 @param url URL to an HDF5 file. File will be created if it doesn't already exist.
 @param dataSet Data set name in HDF5. If the data set already exists, write will fail.
 @param OUT NSError instance if error.
 
 @return YES if successful, NO otherwise.
 */
- (BOOL)writeToHDF5URL:(NSURL*)url dataSet:(NSString*)dataSet error:(NSError**)err;
@end
