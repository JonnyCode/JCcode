//
//  ADMBaseAcquirinoWrapper.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/13/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//


@interface ADMBaseAcquirinoWrapper : NSObject {
@protected
  NSString* filePath;
}

- (id) initWithFile: (NSString *) path;
- (void) throwIfAcquirinoError: (long) acquirinoCode;

@end

