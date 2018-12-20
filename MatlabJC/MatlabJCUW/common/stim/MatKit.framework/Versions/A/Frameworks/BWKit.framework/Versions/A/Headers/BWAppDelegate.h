//
//  BWAppDelegate.h
//  BWKit
//
//  Created by Barry Wark on 11/1/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <BWKit/BWPluginManager.h>

OBJC_EXPORT NSString *BWAppDelegateNewTicketURL;
OBJC_EXPORT NSString *BWAppDelegatePrefsNibName;
OBJC_EXPORT NSString *BWAppDelegateSplashNibName;
OBJC_EXPORT NSString *BWAppDelegateRedirectLogsURLKey;
OBJC_EXPORT NSString *BWAppDelegateShowDebugLogKey;
OBJC_EXPORT NSString *BWAppDelegateAppShouldOpenUntitledDocumentKey;

@interface BWAppDelegate : NSObject {
    IBOutlet NSWindow *splashWindow;
    IBOutlet NSPanel *prefsPanel;
    IBOutlet NSTreeController *pluginsController;
    IBOutlet NSPanel *pluginsPanel;
    IBOutlet NSProgressIndicator *splashProgress;
    
    NSWindowController *splashController;
    
    NSOperationQueue *operationQueue;
    
    NSURL *newTicketURL;
    
    BOOL singleOpenDocument;
    NSString *splashStatus;
}

@property (retain,readwrite) NSURL *newTicketURL;
@property (retain,readonly) NSOperationQueue *operationQueue;
@property (readonly) BWPluginManager *pluginManager;
@property (retain,readonly) NSString *splashStatus;
@property (nonatomic,readonly) NSUInteger buildNumber;

/*!
 @abstract Policy: single open document (in document based apps)
 */
@property (assign,readwrite) BOOL singleOpenDocument;

/*!
 @method showPreferences:(id)sender
 @abstract Present the app prefs. dialog, updating the user defaults only NSOKButton is used to dismiss the dialog.
 @discussion Loads and presents from BWAppDelegatePrefsNibName with self as owner. Make sure to connect self.prefsPanel IBOutlet.
 */
- (IBAction)showPreferences:(id)sender;
- (IBAction)preferencesOK:(id)sender;
- (IBAction)preferencesCancel:(id)sender;

/*!
 @method     
 @abstract   File a new (Trac) ticket.
 @discussion Opens the URL in Info.plist[BWAppDelegateNewTicketURL] in the system browser.
 */
- (IBAction)fileNewTicket:(id)sender;


/*!
 @method     
 @abstract   Load a new user-selected plugin
 @discussion Displays an OpenPanel for selection of plugin.
 */
- (IBAction)loadNewPlugin:(id)sender;

/*!
 @method     
 @abstract   Registered callback for plugin load notifications.
 @discussion Sub-classes should override postLoadPlugin: to perform any necessary processing on loaded plugins.
 */
- (void)postLoadPlugin:(NSNotification*)notificaiton;

/*!
    @method     
    @abstract   Application delegate method called before application terminates
    @discussion Shuts down logging system
*/

- (void)applicationWillTerminate:(NSNotification *)aNotification;

/*!
 @method     
 @abstract   Display splash screen and load plugins.
 @discussion Subclasses should do their own processing, then call super's 
 implementation. Subclasses should override postLoadPlugin: to detect plugin loads. 
 Initializes logging system. Calls initPreferences.
 @seealso postLoadPlugin:
 */
- (void)applicationWillFinishLaunching:(NSNotification*)aNotification;

/*!
 @method     
 @abstract   Remove splash screen.
 */
- (void)applicationDidFinishLaunching:(NSNotification*)aNotification;

/*!
 @method     
 @abstract   Allowed extensions for application plugins.
 @discussion Sub-classes may override to append their own plugin extensions. Default is {@"plugin"}.
 */
- (NSSet*)pluginExtensions;

// !!!:barry:20080519 TODO
- (void)initPreferences;

/*!
 @method     
 @abstract   Dictionary of preferences defaults.
 @discussion If DefaultPrefValues.plist is present in the app's bundle, this
 method loads and returns that plist. Otherwise returns an empty dictionary.
 @seealso initPreferences
*/
- (NSDictionary*)preferencesDefaults;

/**
 Application build number.
 
 Used by crash reporter, and available for rest of app.
 */
- (NSUInteger)buildNumber;
@end

/*
 TODO
 - show PDF/URL
 */
