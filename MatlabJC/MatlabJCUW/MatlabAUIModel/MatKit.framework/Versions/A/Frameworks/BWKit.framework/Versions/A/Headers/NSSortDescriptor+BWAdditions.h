//
//  NSSortDescriptor+BWAdditions.h
//  VizUI
//
//  Created by Barry Wark on 7/16/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSSortDescriptor (BWAdditions)
+ (NSSortDescriptor*)sortDescriptorWithKey:(NSString*)key ascending:(BOOL)ascending;
@end
