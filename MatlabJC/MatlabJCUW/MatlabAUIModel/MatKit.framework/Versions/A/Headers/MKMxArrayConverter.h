//
//  MKMxArrayConverter.h
//  MatKit
//
//  Created by Daniel Carleton on 6/10/08.
//  Copyright 2008 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "mat.h"

// any struct mxArray having this field will be expected to also have an
// timeIntervalSince1970 field that contains a value convertable to an NSDate
// via dateWithTimeIntervalSince1970.  this (kludgy) hint is required because
// object instances can't be passed into mex land.
extern NSString* const MKEpochTimeStructSigField;

@interface MKMxArrayConverter : NSObject {
}

+ (id) objectForArray: (mxArray *) array;

@end

