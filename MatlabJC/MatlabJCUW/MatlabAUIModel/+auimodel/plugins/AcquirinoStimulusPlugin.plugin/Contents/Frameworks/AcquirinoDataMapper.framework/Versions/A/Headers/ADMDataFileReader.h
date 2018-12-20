//
//  ADMDataFileReader.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/10/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import "ADMBaseAcquirinoWrapper.h"
#import "ADMEpochEnumerator.h"

@interface ADMDataFileReader : ADMBaseAcquirinoWrapper {
  ADMItcDevice itcDevice;
}

+ (id) readerForFile: (NSString *) path
         usingDevice: (ADMItcDevice) itcDevice;
- (NSString *) getComment;
- (NSNumber *) getVersionNumber;
- (ADMEpochEnumerator *) epochEnumerator;
- (id) initWithFile: (NSString *) path
        usingDevice: (ADMItcDevice) itcDevice;

@end

