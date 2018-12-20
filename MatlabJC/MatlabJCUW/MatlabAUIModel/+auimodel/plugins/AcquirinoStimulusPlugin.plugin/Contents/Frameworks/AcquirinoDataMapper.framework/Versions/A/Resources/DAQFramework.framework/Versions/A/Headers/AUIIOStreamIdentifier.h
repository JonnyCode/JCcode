//
//  AUIIOStreamIdentifier.h
//  AcqUI
//
//  Created by Barry Wark on 10/19/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIIOStream.h>

/*
 AUIIOStreamIdentifier is an archiveable unique identifier of a stream.
 */
@interface AUIIOStreamIdentifier : NSObject <NSCoding> {
    StreamType type;
	unsigned long channelNumber;
    id controllerID;
}

@property (retain) id controllerID;
@property (assign) unsigned long channelNumber;
@property (assign) StreamType type;

+ (AUIIOStreamIdentifier*)streamIdentifierWithChannelNumber:(unsigned long)channel
                                                       type:(StreamType)streamType
                                               controllerID:(id)aControllerID;

- (id)initWithChannelNumber:(unsigned long)channel
                       type:(StreamType)streamType
               controllerID:(id)aControllerID;

/*!
 @abstract Find the current stream corresponding to this contain'er identifier.
 @discussion Return the true stream corresponding to this container's identifier by searching the registered controller,
 then the DAQ framework's current controller (if it exists). If the controller isn't available, or the stream
 does not exist, return an unbound stream.
 @result AUIIOStream belonging to the appropriate controller, or an unbound stream if no appropriate controller exists.
 */
- (id)stream;

@end
