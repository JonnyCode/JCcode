//
//  NSString+TempFile.h
//  BWKit
//
//  Created by Barry Wark on 4/7/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSString (TempFile)
/*!
 @method     
 @abstract   Build a new temporary file path (in NSTemporaryDirectory())
 @discussion Uses tempnam to create a new file path in NSTemporaryDirectory(). 
 The resulting path is unique at the time of its creation. There is a possible
 race condition, however: an attacker could attempt to create a file at the
 same path between when the path is determined and when the file is opened. Use
 tmpfile or mkstemp when possible.
 @returns NSString containing a an unused path in NSTemporaryDirectory().
*/

+ (NSString*)tempFilePath;
@end
