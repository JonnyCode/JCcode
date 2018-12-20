//
//  NSObject+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 7/31/07.
//  Copyright 2007 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSObject (BWAdditions)
- (void)fireChangeNotificationForKey:(NSString*)key;
@end
