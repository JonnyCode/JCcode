//
//  NSInvocation+BWAdditions.h
//  BWKit
//
//  Created by Barry Wark on 11/29/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface NSInvocation (BWAdditions)
/*!
    @method     
    @abstract   Create an invocation with empty argument list.
    @discussion setArguments: must be called to set the arguments if sel requires arguments.
*/
+ (NSInvocation*)invocationWithTarget:(id)invocationTarget
                             selector:(SEL)sel;

/*!
 @method     
 @abstract   Create an invocation with argument list.
 @discussion All arguments must be type id.
 @param arg_list If nil, arguments must be filled in later.
*/

+ (NSInvocation*)invocationWithTarget:(id)invocationTarget
                 selectorAndArguments:(SEL)sel,...;


/*!
    @method     
    @abstract   Convenience method for calling from other methods
    @discussion <#(comprehensive description)#>
*/

+ (NSInvocation*)invocationWithTarget:(id)invocationTarget
                             selector:(SEL)sel
                            arguments:(va_list)arg_list;
/*!
 @method     
 @abstract   Set arguments for this invocation.
 @discussion Arguments must all be type id.
 @param firstArg Index (in self.methodSignature) of the first argument in arg_list. Index is from 2 (arguments 0 and 1 are self and _cmd). A method that takes one explicit arguemtn (eg. mySelector:) should have firstArg == 0.
*/

- (void)setArguments:(id)arg1,...;
- (void)setArgumentList:(va_list)arg_list
               firstArg:(NSUInteger)firstArg;
@end
