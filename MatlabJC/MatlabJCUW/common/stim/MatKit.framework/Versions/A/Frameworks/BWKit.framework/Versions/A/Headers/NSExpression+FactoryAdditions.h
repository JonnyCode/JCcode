//
//  NSExpression+FactoryAdditions.h
//  BWKit
//
//  Created by Barry Wark on 10/23/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSExpression (FactoryAdditions)
+ (NSExpression*)valueForKeyPathFunctionExpressionForOperand:(NSExpression*)operand
                                                     keyPath:(NSString*)keyPath;
+ (NSExpression*)valueForKeyPathFunctionExpressionForOperand:(NSExpression*)operand
                                               constantValue:(NSString*)arg;
@end
