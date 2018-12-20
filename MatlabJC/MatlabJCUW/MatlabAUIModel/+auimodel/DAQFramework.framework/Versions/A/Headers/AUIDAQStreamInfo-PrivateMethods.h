/*
 *  AUIDAQStreamInfo-PrivateMethods.h
 *  AcqUI
 *
 *  Created by Barry Wark on 1/23/07.
 *  Copyright 2007 Barry Wark. All rights reserved.
 *
 */

#import "DAQFramework/AUIDAQStream.h"

@interface AUIDAQStream(PrivateMethods)
- (void)saveLastSample;
- (NSMutableData *)data;
- (void)setData:(NSMutableData *)newData;
@end

