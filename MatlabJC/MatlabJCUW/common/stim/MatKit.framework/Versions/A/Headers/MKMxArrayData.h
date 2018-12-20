//
//  MKMxArrayData.h
//  MatKit
//
//  Created by Daniel Carleton on 1/23/08.
//  Copyright 2008 __MyCompanyName__. All rights reserved.
//

#import <CoreFoundation/CoreFoundation.h>

extern NSString* const MKMatlabErrorDomain;

/*!
 * @class MKMxArrayData
 * @abstract Concrete NSData subclass for managing mxArray data.
 * @discussion All copy and freeWhenDone parameters are ignored.  Data is not
 * copied, and is automatically freed with mxDestroyArray on dealloc.
 */
@interface MKMxArrayData : NSData {
  NSData* embeddedData;
  bool destroyOnRelease;
}

@property (readonly) NSData* embeddedData;
@property (assign) bool destroyOnRelease;

+ (MKMxArrayData *) dataWithMXBytes: (void *) array;
- (id) initWithMXBytes: (void *) array;
- (id) initWithMXBytes: (void *) array;

/*!
 * @method writeToFile:withName:error:
 * @abstract Writes this Matlab array to a mat file.
 */
- (BOOL) writeToFile: (NSString *) path
      asVariableName: (NSString *) name
               error: (NSError **) error;
@end

