//
//  AUIExternalDeviceManager.h
//  AcqUI
//
//  Created by Barry Wark on 9/28/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@protocol AUIExternalDevice;

@interface AUIExternalDeviceManager : NSObject {
    NSMutableSet *registeredSet;
}

+ (NSSet*)keyPathsForValuesAffectingRegisteredDeviceClasses;
+ (NSSet*)keyPathsForValuesAffectingRegisteredDeviceClassNames;

+ (AUIExternalDeviceManager*)sharedManager;
+ (id<AUIExternalDevice>)nullDevice;

- (void)registerDeviceClass:(Class)d;
- (NSSet*)registeredDeviceClasses;
- (NSSet*)registeredDeviceClassNames;
@end
