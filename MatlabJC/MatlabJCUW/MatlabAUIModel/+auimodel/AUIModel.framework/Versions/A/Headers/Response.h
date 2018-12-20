//
//  Response.h
//  AcqUI
//
//  Created by Barry Wark on 12/29/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import "IOBase.h"

@class Epoch;
@class BWNumericData;

@interface Response :  IOBase  
{
}

- (BWNumericData *)data;
- (void)setData:(NSData *)value;

- (Epoch *)epoch;
- (void)setEpoch:(Epoch *)value;

- (NSManagedObject *)dataContainer;
- (void)setDataContainer:(NSManagedObject *)value;

@end
