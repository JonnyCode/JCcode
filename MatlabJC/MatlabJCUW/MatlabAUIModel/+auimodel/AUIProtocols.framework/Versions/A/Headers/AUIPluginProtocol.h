/*
 *  AUIPluginProtocol.h
 *  AcqUI
 *
 *  Created by Barry Wark on 11/11/06.
 *  Copyright 2006 Barry Wark. All rights reserved.
 *
 */

#include <Cocoa/Cocoa.h>
#include <BWKit/BWKit.h>

#define AUI_STIMULUS_PLUGIN_TYPE @"edu.washington.bwark.acqui.stimulus"
#define AUI_PROTOCOL_PLUGIN_TYPE  @"edu.washingon.bwark.acqui.protocol"
#define AUI_DAQ_CONTROLLER_PLUGIN_TYPE @"edu.washington.bwark.acqui.DAQController"
#define AUI_EXTERNAL_DEVICES_PLUGIN_TYPE @"edu.washington.bwark.acqui.ExternalDevice"

#define AUI_PLUGIN_TYPE_KEY @"AUIPluginType"

@protocol AUIPlugin <BWPlugin>
@end