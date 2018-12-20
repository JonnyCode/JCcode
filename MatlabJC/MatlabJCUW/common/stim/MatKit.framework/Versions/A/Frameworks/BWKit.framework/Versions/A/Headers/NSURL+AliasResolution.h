//
//  NSURL+AliasResolution.h
//  BWKit
//
//  Created by Barry Wark on 10/31/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSURL (BWAliasResolution)
/**
 Creates a URL to a fileURL by resolving all aliases in the path.
 
 CFURLGetFSRef cannot get FSRefs for URLs containing aliases. This method
 iterates the path components, resolving them as needed.
 */
- (NSURL*)urlByResolvingAliasesInPath;
@end
