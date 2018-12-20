//
//  QTMovie+CIImageFrameEnumeration.h
//  BWKit
//
//  Created by Barry Wark on 4/3/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <QTKit/QTKit.h>

@interface QTMovie (CIImageFrameEnumeration)
/*!
 @method     
 @abstract   Return an enumerator of CIImages for each frame of movie.
 @discussion If sampleSize is not a multiple of movie duration, enumerates
 floor(movieDuration/sampleSize) frames.
 @param sampleRate Sample rate (Hz). Must be >= 0.
*/

- (NSEnumerator*)CIImageFrameEnumeratorForSampleRate:(double)sampleRate;
@end
