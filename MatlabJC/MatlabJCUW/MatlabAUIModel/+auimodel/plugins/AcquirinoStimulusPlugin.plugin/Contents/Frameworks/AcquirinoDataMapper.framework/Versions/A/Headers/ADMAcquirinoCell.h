//
//  ADMAcquirinoCell.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 12/20/07.
//  Copyright 2007 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface ADMAcquirinoCell : NSObject {
@private
  NSDate* startDate;
  NSString* dataFileBaseName;
  NSArray* epochFamilies;
  NSString* label;
}

@property (readwrite,retain) NSDate* startDate;
@property (readwrite,retain) NSString* dataFileBaseName;
@property (readwrite,retain) NSArray* epochFamilies;
@property (readwrite,retain) NSString* label;

@end

