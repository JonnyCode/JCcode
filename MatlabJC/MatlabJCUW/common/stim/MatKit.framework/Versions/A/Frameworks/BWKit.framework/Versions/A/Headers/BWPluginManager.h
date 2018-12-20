//
//  BWPluginManager.h
//  VizUI
//
//  Created by Barry Wark on 10/13/06.
//  Copyright 2006 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>

OBJC_EXPORT NSString *BWPluginLoadedNotification;
OBJC_EXPORT NSString *BWPluginUnloadedNotificaiton;

// keys in NSNotification.userInfo for the BWPlugin*Notification
OBJC_EXPORT NSString *BWPluginClassKey;
OBJC_EXPORT NSString *BWPluginTypeKey;
OBJC_EXPORT NSString *BWPluginClassNameKey;

@interface BWPluginManager : NSObject {
    NSMutableSet *pluginExtensions;
    NSMutableDictionary *plugins; //pluginType => {set of plugins}
    
    NSString *status;
    __weak NSTextField *statusField; //will be refreshed at each status change if not nil. weak ref.
}

@property(retain) NSString *status;
@property(assign) __weak NSTextField *statusField;
@property(retain) NSMutableSet *pluginExtensions;
@property(retain,readonly) NSMutableDictionary *plugins;

+ (BWPluginManager*)manager;
- (void)loadPluginsAtPath:(NSString*)path;
- (void)registerPluginClass:(Class)pluginClass ofType:(NSString*)type;
- (void)loadedNewPluginClass:(Class)cls ofType:(NSString*)type; //called for each loaded plugin principal class. Posts BWPluginLoadedNotification, with class and type of plugin in userInfo dictionary. subclasses should register for particular types if they want to do post-load processing.


- (void)clearLoadedPlugins;

- (NSSet*)pluginsOfType:(id)type includeDecprecated:(BOOL)includeDecprecated;
- (NSSet*)pluginsOfType:(id)type;

- (NSSet*)pluginsForProtocol:(Protocol*)p includeDecprecated:(BOOL)includeDecprecated;
- (NSSet*)pluginsForProtocol:(Protocol*)p;
@end
