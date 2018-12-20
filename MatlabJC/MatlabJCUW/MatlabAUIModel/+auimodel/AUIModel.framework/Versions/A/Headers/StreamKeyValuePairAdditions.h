//
//  StreamKeyValuePairAdditions.h
//  AcqUI
//
//  Created by Barry Wark on 1/1/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <AUIModel/StreamKeyValuePair.h>

@interface StreamKeyValuePair (Additions)
/*!
 @method     
 @abstract   Convert self.streamValue to an AUIIOStream
 @result AUIIOStream (unbound if AUIIOStreamIdentifier.stream can't find a current controller.)
 */

+ (id)kvpValueForInstance:(NSManagedObject*)kvp;
@end
