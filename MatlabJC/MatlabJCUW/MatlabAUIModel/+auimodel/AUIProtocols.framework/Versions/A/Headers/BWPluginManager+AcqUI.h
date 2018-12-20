//
//  BWPluginManager+AcqUI.h
//  AcqUI
//
//  Created by Barry Wark on 7/16/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWKit.h>
#import <AUIProtocols/AUIStimulusGeneratorProtocol.h>
#import <AUIProtocols/AUIStimulusProtocol.h>
#import <AUIProtocols/AUIPluginProtocol.h>

@interface BWPluginManager (AcqUI)
- (Class)stimulusClassWithID:(NSString*)stimulusID;
- (id<AUIStimulusProtocol>)protocolWithID:(NSString*)protocolID;
@end
