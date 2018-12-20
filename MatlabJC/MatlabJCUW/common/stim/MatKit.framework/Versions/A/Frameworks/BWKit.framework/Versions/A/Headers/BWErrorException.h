//
//  BWErrorException.h
//  BWKit
//
//  Created by Barry Wark on 11/17/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>

OBJC_EXPORT NSString * const BWErrorExceptionName;
OBJC_EXPORT NSString * const BWErrorExceptionErrorKey;

/*!
    @class
    @abstract    NSException subclass that wraps an NSError
    @discussion  Stores the NSError instance in error property. userInfo is error's userInfo.
*/

@interface BWErrorException : NSException {
    NSError *error;
}

@property (copy,readonly) NSError *error;

+ (BWErrorException*)exceptionWithError:(NSError*)err;
+ (void)raiseError:(NSError*)err;
- (id)initWithError:(NSError*)error;
@end
