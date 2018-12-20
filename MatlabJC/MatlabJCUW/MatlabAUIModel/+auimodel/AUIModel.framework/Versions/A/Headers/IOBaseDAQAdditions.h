//
//  IOBaseDAQAdditions.h
//  VizUI
//
//  Created by Barry Wark on 10/13/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <AUIModel/IOBase.h>
#import <DAQFramework/AUIExternalDevice.h>
#import <BWKit/BWNumericDataType.h>

@class AUIIOStream;

@interface IOBase (IOBaseDAQAdditions)
- (NSDictionary*)streamProperties;
- (id<AUIExternalDevice>)externalDevice;
- (double)samplingRate;

- (NSString*)streamInfoUserDescription;
- (NSString*)externalDeviceUnits;

- (NSString*)dtypeString; //numpy-style dtype string for this IOBases' data
+ (NSString*)dtypeStringForInstance:(IOBase*)ioBase;
- (BWNumericDataType*)dtype;
+ (BWNumericDataType*)dtypeForInstance:(IOBase*)ioBase;
@end
