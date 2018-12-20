//
//  NumberKeyValuePair.h
//  AcqUI
//
//  Created by Barry Wark on 12/19/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import <BWKit/KeyValuePair.h>


@interface NumberKeyValuePair :  KeyValuePair  <KeyValuePairSubEntityMethods>
{
}

+ (NSEntityDescription*)entitySubclassForValue:(NSValue*)n
                        inManagedObjectContext:(NSManagedObjectContext*)moc;
@end


