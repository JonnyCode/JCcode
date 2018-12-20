//
//  ADMEpochEnumerator.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/11/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ADMBaseAcquirinoWrapper.h"
#import "ADMAcquirinoEpoch.h"

typedef enum {
  ADMItc16Device = 0,
  ADMItc18Device
} ADMItcDevice;

@interface ADMEpochEnumerator : ADMBaseAcquirinoWrapper {
@private
  ADMAcquirinoEpoch* lastEpoch;
  long epochCount;
  long currentEpoch;
  long maxPointsPerEpoch;
  BWNumericDataType* gclampWaveformDtype;
  ADMItcDevice itcDevice;
}

- (id) initWithEpochCount: (long) count
                  forFile: (NSString *) path
              usingDevice: (ADMItcDevice) device;

- (ADMAcquirinoEpoch *) nextObject;
- (NSArray *) allObjects;

@end

