/*
 *  BWPluginProtocol.h
 *  VizUI
 *
 *  Created by Barry Wark on 7/16/07.
 *  Copyright 2007 Barry Wark. All rights reserved.
 *
 */

#include <Cocoa/Cocoa.h>

@protocol BWPlugin <NSObject>
+(id)pluginType;
@end

