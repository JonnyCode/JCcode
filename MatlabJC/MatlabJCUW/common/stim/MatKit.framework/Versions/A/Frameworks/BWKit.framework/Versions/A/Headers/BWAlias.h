//
//  BWAlias.h
//  BWKit
//
//  Created by Barry Wark on 10/31/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BDAlias.h>

@interface BWAlias : BDAlias {

}

+ (BWAlias*)aliasWithPath:(NSString *)path relativeToPath:(NSString *)relPath;

- (NSURL*)URLRelativeToURL:(NSURL*)root;
@end
