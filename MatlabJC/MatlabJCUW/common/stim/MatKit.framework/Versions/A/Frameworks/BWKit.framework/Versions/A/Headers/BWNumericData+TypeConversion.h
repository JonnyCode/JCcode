//
//  BWNumericData+TypeConversion.h
//  BWKit
//
//  Created by Barry Wark on 4/30/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWNumericData.h>

@interface BWNumericData (TypeConversion)
- (BWNumericData*)dataByConvertingToType:(BWDataType)newDataType
                             sampleBytes:(NSUInteger)newSampleBytes
                               byteOrder:(CFByteOrder)newByteOrder;
@end
