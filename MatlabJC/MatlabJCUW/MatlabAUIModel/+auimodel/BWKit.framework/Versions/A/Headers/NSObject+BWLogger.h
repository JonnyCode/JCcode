//
//  NSObject+BWLogger.h
//  BWKit
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWLogger.h>


@interface NSObject (BWLoggerAdditions)

- (BWLogger*)logger;
+ (BWLogger*)logger;

@end