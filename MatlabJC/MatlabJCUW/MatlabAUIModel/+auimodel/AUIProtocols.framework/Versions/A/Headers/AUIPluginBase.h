//
//  AUIPluginBase.h
//  AcqUI
//
//  Created by Barry Wark on 1/18/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "AUIProtocols/AUIPluginProtocol.h"
#import <BWKit/BWKit.h>


@interface AUIPluginBase : BWPluginBase <AUIPlugin> {

}

+(id)pluginType;
+(id)pluginID;
@end
