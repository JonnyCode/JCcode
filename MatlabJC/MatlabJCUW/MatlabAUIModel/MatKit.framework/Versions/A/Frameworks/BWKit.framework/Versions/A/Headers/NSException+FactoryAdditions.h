//
//  NSException+FactoryAdditions.h
//  BWKit
//
//  Created by Barry Wark on 11/17/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSException (FactoryAdditions)
+ (void)raiseGenericFormat:(NSString*)fmt,...;
+ (void)raiseNotImplemented;
@end
