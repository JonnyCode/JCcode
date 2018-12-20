//
//  ADMAbstractNotebookFileParser.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/21/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

// bring in the model classes our client may need
#import "ADMAcquirinoCell.h"
#import "ADMAcquirinoEpochFamily.h"
#import "ADMAcquirinoNote.h"

@interface ADMAbstractNotebookFileParser : NSObject {
  NSString* filePath;
}

// getters for these properties are defined using the key-value coding support
// in pyobjc
@property (readonly,copy) NSString* filePath;
@property (readonly,copy) NSString* experimentPurpose;
@property (readonly,copy) NSDate* experimentStartDate;
@property (readonly,assign) NSArray* notes;
@property (readonly,assign) NSArray* cells;

+ (id) parserForFile: (NSString *) path;
- (id) initWithFile: (NSString *) path;
- (void) parse;

@end

