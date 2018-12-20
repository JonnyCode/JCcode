//
//  BWNumericData.h
//  BWKit
//
//  Created by Barry Wark on 4/8/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWNumericDataType.h>

OBJC_EXPORT NSString *BWNumericDataException;

@interface BWNumericData : NSData {
    NSData *data;
    BWNumericDataType *dtype;
    NSArray *shape; //array of dimension shapes (NSNumber<unsigned>)
}

@property (retain,readonly) BWNumericDataType* dtype;
@property (copy,readonly) NSArray* shape;
@property (readonly) NSUInteger ndims;
@property (readonly) NSUInteger nSamples; //number of samples of dtype
@property (readonly) BWDataType dataType;
@property (readonly) NSUInteger sampleBytes;
@property (readonly) CFByteOrder byteOrder;


+ (BWNumericData*)numericDataWithData:(NSData*)theData
                                dtype:(BWNumericDataType*)_dtype
                                shape:(NSArray*)shapeArray;

/*!
 @method     
 @abstract   DESIGNATD INITIALIZER. Initialize a BWNumericData object from data, dtype, and shape
 @discussion A BWNumericData instance can be initialized from a numpy array (in python)::
 numericData = BWNumericData.alloc().initWithData_dtypeString_shape_(numpy_array,
 										numpy_array.dtype.str,
 										numpy_array.shape)
 
 @throws NSException if shape is non-nil and the product of the shape elements does not match
 the data size (length/sampleBytes(dtype)).
 @param theData NSData* (copied)
 @param dtype dtype (retained)
 @param shapeArray data shape (ala numpy) (copied). may be nil for 1-D data.
*/

- (id)initWithData:(NSData*)theData
             dtype:(BWNumericDataType*)_dtype
             shape:(NSArray*)shapeArray;

- (id)initWithData:(NSData*)theData
             dtypeString:(NSString*)dtypeString
             shape:(NSArray*)shapeArray;

- (NSUInteger)nSamples;
- (BWDataType)dataType;
- (NSUInteger)sampleBytes;
- (CFByteOrder)byteOrder;


- (void*)samplePointer:(NSUInteger)sample;
- (NSNumber*)sampleValue:(NSUInteger)sample;
@end
