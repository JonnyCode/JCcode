//
//  AUIIOStreamIdentifierContainer.h
//  AcqUI
//
//  Created by Barry Wark on 10/20/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIIOStream.h>

@class AUIIOStreamIdentifier;

@interface AUIIOStreamIdentifierContainer : NSObject <NSCoding> {
    id streamIdentifier;
}

@property (retain) id streamIdentifier;

- (id)initWithStreamIdentifier:(AUIIOStreamIdentifier*)identifier;
- (id)initWithStream:(AUIIOStream*)stream;

- (unsigned long)channelNumber;
- (StreamType)type;
- (id)controllerID;

- (id)stream;

@end
