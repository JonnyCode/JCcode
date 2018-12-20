/*
 *  BWCocoaMPLPlotViewBindings.h
 *  AcqUI
 *
 *  Created by Barry Wark on 6/18/07.
 *  Copyright 2007 Barry Wark. All rights reserved.
 *
 */

@protocol BWCocoaMPLPlotViewBindings <NSObject>
/*
 Conforming classes must be KVO compliant for these keys:
     lengthArray : {array of float}
     multiAxes : {BOOL}
     unitsArray : {array of strings}
     xlabel : {string}
     ylabelArray : {array of strings}
     dataAndTypeArrays : {array[2]}
         array[0] is a list of data instances
         array[1] is a list of numpy.dtype (or anything that can be used to construct a dtype, e.g. a numpy.dtype.str)
 
 All arrays much have same length, and common indexing (i.e. dataAndTypeArrays[0] has units unitsArray[0] etc.)
*/

@property (retain,readonly) NSArray *lengthArray;
@property (assign,readwrite) BOOL multiAxes;
@property (retain,readonly) NSArray *unitsArray;
@property (retain,readwrite) NSString *xlabel;
@property (retain,readwrite) NSArray *ylabelArray;
@property (retain,readonly) NSArray *dataAndTypeArrays;
@end