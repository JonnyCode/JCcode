//
//  NSDictionary+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 7/22/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSDictionary (BWAdditions)
- (BOOL)hasKey:(id)key;
- (BOOL)hasValuesForKeys:(NSArray*)keys;
- (NSInteger)intForKey:(id)key;
- (BOOL)boolForKey:(id)key;
- (id)objectForKey:(id)key defaultValue:(id)defaultValue; //returns value for key or default value if key not in self.
@end
