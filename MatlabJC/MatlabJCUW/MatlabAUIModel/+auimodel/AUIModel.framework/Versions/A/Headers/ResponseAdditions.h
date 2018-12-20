//
//  ResponseAdditions.h
//  AcqUI
//
//  Created by Barry Wark on 12/29/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <AUIModel/Response.h>
#import <AUIModel/CachedData.h>
#import <BWKit/BWKit.h>

@interface Response (ResponseAdditions) <CachedData>
- (BWNumericData*)data;
- (void)setData:(NSData*)value;
- (NSTimeInterval)duration;

- (void)flushCachedData;
- (void)saveResponseDataToHDF5;
+ (void)saveResponseData:(NSManagedObject*)response toResponseDataFileURL:(NSURL*)url;
@end
