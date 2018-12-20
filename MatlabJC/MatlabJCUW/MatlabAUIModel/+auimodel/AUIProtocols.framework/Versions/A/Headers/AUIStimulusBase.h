//
//  AUIStimulusBase.h
//  AcqUI
//
//  Created by Barry Wark on 2/26/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "AUIStimulusGeneratorProtocol.h"
#import <AUIProtocols/AUIPluginProtocol.h>

/*!
AUIStimulusBase provides implementations or stubs for all
 required methods of AUIStimulus.
 
 Stimulus generators should probably subclass from this
 class.
 */
@interface AUIStimulusBase : NSObject <AUIStimulus, AUIPlugin> {
	double sampleRate;
	NSTimeInterval length;
}

@property (assign,readwrite) NSTimeInterval length;
@property (assign,readwrite) double sampleRate;

/*!
Implemented as [self stimulusDataForParamsDictionary:[self paramsDictionary]];
 */
-(NSData*)stimulusData;

/*!
Calls [self stimulusDataForParamsDictionary:p version:[self version]];
 */

-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p;

/*!
 Stub implementation (will throw an exception if called).
 */
-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p version:(NSUInteger)version;

/*!
 For backwards compatibility, calls [self stimulusDataForParamsDictionary:version:] (ignoring fileSystemResources and subStimuli).
 Future subclasses should override this method.
 */
-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p version:(NSUInteger)version fileSystemResources:(NSSet*)fileSystemResources subStimuli:(NSSet*)subStimuli;

/*!
Stub implementation (will throw an exception if called).
 */
-(NSDictionary*)paramsDictionary;

/*!
Stub implementation (will throw an exception if called).
 */
-(NSUInteger)version;

/*!
Returns the Class' bundle's ID. This is probably good enough for a stimulusID.
 */
+(NSString*)stimulusID;

+(id)pluginType;
@end
