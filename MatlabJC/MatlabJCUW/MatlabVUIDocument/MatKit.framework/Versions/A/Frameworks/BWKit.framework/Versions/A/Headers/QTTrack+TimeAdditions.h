//
//  QTTrack+TimeAdditions.h
//  BWKit
//
//  Created by Barry Wark on 4/3/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <QTKit/QTKit.h>


@interface QTTrack (TimeAdditions)
- (long long)sampleTime;
- (long)timeScale;
- (long long)numSamples;
@end
