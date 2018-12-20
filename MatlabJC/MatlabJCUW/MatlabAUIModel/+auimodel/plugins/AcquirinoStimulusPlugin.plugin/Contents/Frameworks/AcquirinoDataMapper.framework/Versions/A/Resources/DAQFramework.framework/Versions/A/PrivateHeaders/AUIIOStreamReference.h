//
//  AUIIOStreamReference.h
//  AcqUI
//
//  Created by Barry Wark on 10/19/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIIOStream.h>
#import <DAQFramework/AUIIOStreamIdentifierContainer.h>

/*
 AUIIOStreamReference encapsulates the unique identifying information for an AUIIOStreamInfo. It conforms to the 
 NSCoding protocol and can be encoded directly. On decoding, however, it will decode to the stream info it
 references (i.e. the AUIIOStreamInfo corresponding to the same channelNumber and type in the registered
             controller of that controllerID)..
 If there is no registered controller with the given ID, returns the corresponding (channelNumber and type)
 stream from any registered controller that can parse the given controllerID's configData.
 
 If no such channel exists, an AUIIOStreamReference decodes to nil.
 */
@interface AUIIOStreamReference : AUIIOStreamIdentifierContainer {
    
}

+ (id)streamReferenceForStream:(AUIIOStream*)stream;

- (id)awakeAfterUsingCoder:(NSCoder *)coder;
@end
