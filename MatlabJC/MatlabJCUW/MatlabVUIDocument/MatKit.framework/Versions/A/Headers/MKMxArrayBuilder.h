//
//  MKMxArrayBuilder.h
//  MatKit
//
//  Created by Daniel Carleton on 1/17/08.
//  Copyright 2008 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <MatKit/MKMxArrayData.h>

extern NSString* const MKMatlabErrorDomain;

/*!
 * @class MKMxArrayBuilder
 * @abstract Encapsulates building Matlab data structures from Foundation types
 * such as NSDictionary, NSArray, NSString, and NSData.
 */
@interface MKMxArrayBuilder : NSObject {

}

/*!
 * @method arrayForObject:
 * @abstract Recursively builds a Matlab array for the given object.
 * @discussion Matlab arrays are built for the following supported types:
 * NSString, NSDictionary, id<NSFastEnumeration>, BWNumericData, and
 * NSManagedObject.  These are transformed into character vectors, structs,
 * matrices, numeric vectors, and structs respectively.
 *
 * Any object is allowed inside of NSFastEnumeration instances, because they're
 * turned into cell arrays.
 *
 * Key-value coding is used to construct structs for NSManagedObjects, which are
 * supported because they're known to conform to the informal NSKeyValueCoding
 * protocol in all cases.  All attributes specified in the attributeRecursion
 * dictionary for NSObject classes are recursed along.
 */
+ (MKMxArrayData *) arrayForObject: (id) object;

/*!
 * @method setAttributeRecursion:
 * @abstract Indicate which attributes to recurse along when buildings structs
 * for NSManagedObjects.
 * @discussion Supply a dictionary of classname => array of attribute names.  If
 * this method isn't called before NSManagedObjects are encountered during
 * array building, an exception is thrown.
 */
+ (void) setAttributeRecursion: (NSMutableDictionary *) classMap;
+ (NSMutableDictionary *) getAttributeRecursion;

@end

