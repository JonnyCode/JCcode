//
//  NSMutableDictionary+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 7/30/07.
//  Copyright 2007 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSMutableDictionary (BWAdditions)
- (id)objectForKey:(id)key settingDefaultValue:(id)defaultValue; //returns value for key, setting value to default if it is not already set.
@end
