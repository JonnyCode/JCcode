//
//  VUIFileSystemResource.h
//  Ovation
//
//  Created by Barry Wark on 10/30/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWFileSystemResource.h>

@class AnalysisRecord;

@interface VUIFileSystemResource : BWFileSystemResource {

}

@property (retain) AnalysisRecord * analysisRecord;

- (NSURL*)relativeRoot;
@end
