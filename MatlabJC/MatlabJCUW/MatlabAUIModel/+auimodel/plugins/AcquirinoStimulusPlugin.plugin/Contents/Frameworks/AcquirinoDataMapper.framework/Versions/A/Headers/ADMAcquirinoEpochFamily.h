//
//  ADMAcquirinoEpochFamily.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/20/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "ADMAcquirinoNote.h"

@interface ADMAcquirinoEpochFamily : NSObject {
  NSRange epochNumberRange;
  NSSet* excludedEpochNumbers;
  NSString* protocolComment;
  ADMAcquirinoNote* tailNote;
  NSDictionary* parameters;
  NSDictionary* epochSettings;
  NSDate* startDate;
}

@property (readwrite,assign) NSRange epochNumberRange;
@property (readwrite,copy) NSString* protocolComment;
@property (readwrite,retain) ADMAcquirinoNote* tailNote;
@property (readwrite,retain) NSDictionary* parameters;
@property (readwrite,retain) NSDictionary* epochSettings;
@property (readwrite,retain) NSDate* startDate;

@end

