//
//  BWAlias.h
//  BWKit
//
//  Created by Barry Wark on 10/31/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/NDAlias.h>

@interface BWAlias : NDAlias {

}

- (NSURL*)URLRelativeToURL:(NSURL*)root;
@end
