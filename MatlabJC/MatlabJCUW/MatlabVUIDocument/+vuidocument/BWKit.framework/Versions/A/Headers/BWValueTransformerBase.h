//
//  BWValueTransformerBase.h
//  BWKit
//
//  Created by Barry Wark on 5/1/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


/*!
 @class
 @abstract    Self-registering NSValueTransformer subclass
 @discussion  Registers self in +load under [self class] name.
*/

@interface BWValueTransformerBase : NSValueTransformer {

}

@end
