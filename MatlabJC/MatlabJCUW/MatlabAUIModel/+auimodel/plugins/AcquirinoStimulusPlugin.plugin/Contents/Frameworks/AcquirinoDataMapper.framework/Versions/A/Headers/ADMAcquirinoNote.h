//
//  ADMAcquirinoNote.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/31/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface ADMAcquirinoNote : NSObject {
  NSString* comment;
  NSDate* timestamp;
}

@property (readwrite,copy) NSString* comment;
@property (readwrite,copy) NSDate* timestamp;

@end

