//
//  NSError+BWAdditions.h
//  BWKit
//
//  Created by Barry Wark on 1/28/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSError (BWAdditions)
/*!
    @method     
    @abstract   Create an NSError object with given failure Reason.
    @discussion Convenience method for NSError creation.
*/

+ (id)errorWithDomain:(NSString *)domain
                 code:(NSInteger)code
      localizedReason:(NSString*)localizedReason;

/*!
    @method     
    @abstract   Create an NSError object with given failure reason and underlying error.
    @discussion Convenience method for NSError creation.
*/

+ (id)errorWithDomain:(NSString *)domain
                 code:(NSInteger)code
      localizedReason:(NSString*)localizedReason
      underlyingError:(NSError*)originalError;


/*!
    @method     
    @abstract   Create an NSError object from an exception e.
    @discussion <#(comprehensive description)#>
*/

+ (id)errorWithException:(NSException*)e;
@end
