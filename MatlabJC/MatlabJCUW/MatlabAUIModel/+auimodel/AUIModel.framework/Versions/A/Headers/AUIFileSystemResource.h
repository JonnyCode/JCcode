//
//  AUIFileSystemResource.h
//  AcqUI
//
//  Created by Barry Wark on 1/8/09.
//  Copyright 2009 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWKit.h>
#import "Experiment.h"

@interface AUIFileSystemResource : BWFileSystemResource {

}

@property (retain) Experiment * experiment;

+ (NSURL*)relativeRootForInstance:(NSManagedObject*)resource;
@end
