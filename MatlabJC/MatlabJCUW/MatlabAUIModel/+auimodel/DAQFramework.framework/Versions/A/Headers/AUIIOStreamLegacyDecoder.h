//
//  AUIIOStreamInfoLegacacyDecoder.h
//  AcqUI
//
//  Created by Barry Wark on 10/22/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIDAQStream.h>

/*
 Old experiment data contains encoded AUIIOStreamInfo. We need to be able to decode the protocol settings
 DAQ settings from these data files. AUIIOStreamInfoLegacyDecoder adds (dumb) decoding for AUIIOStreamInfo
 so that we can decode archives containing AUIIOStreamInfos and get the properties of
 those streams.
 */

@interface AUIIOStreamLegacyDecoder : AUIDAQStream <NSCoding> {

}

@end
