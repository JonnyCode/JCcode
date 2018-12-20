//
//  BWMutableNumericData.h
//  BWKit
//
//  Created by Barry Wark on 4/9/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWNumericDataType.h>

@interface BWMutableNumericData : NSMutableData {
    NSMutableData *data;
    BWNumericDataType *dtype;
    NSArray *shape; //array of dimension shapes (NSNumber<unsigned>)
}

@property (retain,readonly) BWNumericDataType *dtype;
@property (copy,readonly) NSArray* shape;
@property (readonly) NSUInteger ndims;

/*!
 @method     
 @abstract   Initialize a BWNumericData object from data, dtype, and shape.
 @discussion Data retained (not copied)
 @throws NSException if shape is non-nil and the product of the shape elements does not match
 the data size (length/sampleBytes(dtype)).
 */

- (id)initWithData:(NSMutableData*)_data
             dtype:(BWNumericDataType*)_dtype
             shape:(NSArray*)_shape;
@end
