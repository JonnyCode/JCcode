//
//  AUIIOStreamProperties.h
//  AcqUI
//
//  Created by Barry Wark on 10/19/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIIOStream.h>
#import <DAQFramework/AUIIOStreamIdentifier.h>
#import <DAQFramework/AUIIOStreamIdentifierContainer.h>

@interface AUIIOStreamProperties : AUIIOStreamIdentifierContainer {
    NSDictionary *properties;
}

@property (retain) NSDictionary *properties;

+ (id)streamPropertiesForStream:(AUIIOStream*)stream;
+ (id)streamPropertiesForStream:(AUIIOStream*)stream
               streamProperties:(NSDictionary*)theProperties;
+ (id)streamProperitiesForStreamIdentifier:(AUIIOStreamIdentifier*)theIdentifier
                          streamProperties:(NSDictionary*)theProperties;

/*
 Removes "controller", "channelNumber", "type"
 
 @see AUIIOStreamInfo.streamPropertyKeys
 */
/*!
    @method     
    @abstract   Init with a (bound) stream
    @discussion Calculates properties from stream. Unbound streams may not have their properties set. Unbound streams should use initWithStreamIdentifier:streamProperties:
*/

- (id)initWithStream:(AUIIOStream*)stream;

/*!
    @method     
    @abstract   Designated initializer
    @discussion Can be used for unbound streams (to provide a pre-calculatd properties)
*/

- (id)initWithStreamIdentifier:(AUIIOStreamIdentifier*)identifier 
              streamProperties:(NSDictionary*)theProperties;  //designated init

@end
