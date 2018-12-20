//
//  BWOperation.h
//  BWKit
//
//  Created by Barry Wark on 12/24/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface BWOperation : NSOperation {
    NSMutableArray *callbacks; //array of BWOperationCallbackStruct
    
    NSInvocation *invocation;
    
    BOOL callbackOnMainThread;
    
    NSString *name;
    
    NSOperationQueue *queue;
}

@property (retain,readonly) NSMutableArray *callbacks;
@property (retain,readonly) NSInvocation *invocation;
@property (assign,readonly) BOOL callbackOnMainThread; //default = YES
@property (copy,readwrite) NSString *name;
@property (retain,readwrite) NSOperationQueue *queue;

/*!
 @method     
 @abstract   Factory method for constructing a BWOperation
 @discussion God help you if you don't provide the right number of arguments.
 @param selectorAndArguments A SEL that takes arguments of type id.
 @param callbackOnMainThread If YES, all callback/errbacks happen on the main thread.
*/

+ (BWOperation*)operationWithTarget:(id)theTarget
                     callbackTarget:(id)callbackTarget
                   callbackSelector:(SEL)callbackSelector
                      errbackTarget:(id)errbackTarget
                    errbackSelector:(SEL)errbackSelector
               selectorAndArguments:(SEL)theSelector,...;

+ (BWOperation*)operationWithTarget:(id)theTarget
                     callbackTarget:(id)callbackTarget
                   callbackSelector:(SEL)callbackSelector
                      errbackTarget:(id)errbackTarget
                    errbackSelector:(SEL)errbackSelector
               callbackOnMainThread:(BOOL)mainThread
               selectorAndArguments:(SEL)theSelector,...;

/*!
    @method     
    @abstract   Initialize a BWOperation with given targets etc.
    @discussion arg_list must be an initialized va_list. God help you if you don't provide the right number of arguments. arg_list should be va_end'ed after calling this method.
*/

- (id)initWithTarget:(id)theTarget
      callbackTarget:(id)callbackTarget
    callbackSelector:(SEL)callbackSelector
       errbackTarget:(id)errbackTarget
     errbackSelector:(SEL)errbackSelector
callbackOnMainThread:(BOOL)mainThread
            selector:(SEL)theSelector
           arguments:(va_list)arg_list;

/*!
    @method     
    @abstract   Add given callback/errback pair to this operation's callback chain.
    @discussion Attempts to mimick Twisted's callback semantics. // !!!:barry:20080118 document final behavior here
*/

- (void)addCallback:(SEL)callbackSelector
     callbackTarget:(id)callbackTarget
            errback:(SEL)errbackSelector
      errbackTarget:(id)errbackTarget;
/*!
 @method
 @abstract   Operation's main method. Calls this BWOperation's target.selector and callbacks/errbacks on the main thread if self.callbackOnMainThread is True. Otherwise callsback on current thread.
 */

- (void)main;

/*!
    @method     
    @abstract   Cancel this operation
    @discussion Calls self.target.cancel (if implemented), then [super cancel].
*/

- (void)cancel;

/*!
    @method     
    @abstract   Add this operation to an NSOperationQueue
    @discussion Must be used for chained callbacks (i.e. those that return a BWOperation).
*/

- (void)addToQueue:(NSOperationQueue*)theQueue;
@end
