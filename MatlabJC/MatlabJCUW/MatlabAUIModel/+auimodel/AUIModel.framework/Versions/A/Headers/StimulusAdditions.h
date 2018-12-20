//
//  StimulusAdditions.h
//  AcqUI
//
//  Created by Barry Wark on 12/11/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWKit.h>
#import <AUIModel/CachedData.h>
#import "Stimulus.h"

@interface Stimulus (StimulusAdditions) <CachedData>
- (NSDictionary*)parameters;
- (void)setParameters:(NSDictionary*)d;
- (BWNumericData*)reconstructedData;
- (BWNumericData*)reconstructedDataWithPluginManager:(id)pluginManager;
- (void)reconstructData; //sets data to reconstructed data (does not overwrite)

- (void)flushCachedData;
@end
