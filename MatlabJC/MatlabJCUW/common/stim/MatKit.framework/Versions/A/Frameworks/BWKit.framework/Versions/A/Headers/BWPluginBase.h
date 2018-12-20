//
//  BWPluginBase.h
//  VizUI
//
//  Created by Barry Wark on 8/12/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWPluginProtocol.h>

OBJC_EXPORT NSString *BWPluginInfoPlistTypeKey; // Info.plist key for plugin type

@interface BWPluginBase : NSObject <BWPlugin> {

}

+ (id)pluginType;
@end
