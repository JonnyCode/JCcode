//
//  BWNumericDataType.h
//  BWKit
//
//  Created by Barry Wark on 5/3/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>

typedef enum {
    BWUndefinedDataType = 0,
    BWIntegerDataType,
    BWUnsignedIntegerDataType,
    BWFloatingPointDataType,
    BWComplexFloatingPointDataType
} BWDataType;

/*!
    @class
    @abstract    Class to represent the data type of a BWNumericData instance's data.
    @discussion  BWNumericDataType is analogous to NumPy's dtype.
*/

@interface BWNumericDataType : NSObject <NSCoding,NSCopying> {
    BWDataType dataType;
    NSUInteger sampleBytes;
    CFByteOrder byteOrder;
}

@property (assign,readonly) BWDataType dataType;
@property (assign,readonly) NSUInteger sampleBytes;
@property (assign,readonly) CFByteOrder byteOrder;
@property (readonly) NSString *dtypeString;

+ (BWNumericDataType*)dataType:(BWDataType)theType
                   sampleBytes:(NSUInteger)theSampleBytes
                     byteOrder:(CFByteOrder)theByteOrder;

+ (BWNumericDataType*)dataTypeWithDtypeString:(NSString*)dtypeString;

+ (NSString*)dtypeStringForDataType:(BWDataType)dataType
                        sampleBytes:(NSUInteger)sampleBytes
                          byteOrder:(CFByteOrder)byteOrder;

+ (BWDataType)dataTypeForDtypeString:(NSString*)dtypeString;
+ (NSUInteger)sampleBytesForDtypeString:(NSString*)dtypeString;
+ (CFByteOrder)byteOrderForDtypeString:(NSString*)dtypeString;

- (id)initWithDataType:(BWDataType)theType
           sampleBytes:(NSUInteger)theSampleBytes
             byteOrder:(CFByteOrder)theByteOrder;

- (BOOL)isEqualToDataType:(BWNumericDataType*)otherDType;
@end
