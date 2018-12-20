//
//  AUIExternalDeviceBase.h
//  AcqUI
//
//  Created by Barry Wark on 10/5/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "DAQFramework/AUIExternalDevice.h"


@interface AUIExternalDeviceBase : NSObject <AUIExternalDevice, NSCoding> {
    NSString *uuid;
}

- (NSString *)uuid;
- (void)setUuid:(NSString *)anUuid;
@end
