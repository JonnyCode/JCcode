//
//  NSString+RelativePathComputation.h
//  BWKit
//
//  Created by Barry Wark on 11/3/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSString (RelativePathComputation)
- (NSString*)pathRelativeToPath:(NSString*)path;
@end
