//
//  AUIStimulusProtocolBase.h
//  AcqUI
//
//  Created by Barry Wark on 2/26/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "AUIStimulusProtocol.h"
#import "AUIStimulusGeneratorProtocol.h"
#import "EpochDescription.h"
#import "AUIDAQStream.h"
#import "StimulusDescription.h"
#import "ResponseDescription.h"
#import "AUIProtocols/AUIPluginProtocol.h"
#import <AUIProtocols/AUIOnlineAnalystManager.h>

/**
Class that provides implementations and or stubs for the AUIStimulusProtocol, plus
 a few convenience methods to make filling epochs easier.
 
 Protocols should probably subclass from this class.
 
 Subclasses should provide the KVO key name for their settings panel, and their
 settings panel nib names in their Info.plist files with keys:
 AUIProtocolSettingsWindowKVOKey
 AUIProtocolSettingsWindowNibName
 
 Overrides setNilValueForKey: to set 0 for key.
 */
@interface AUIStimulusProtocolBase : NSObject <AUIStimulusProtocol, AUIPlugin>
{
	BOOL isDefault;
    BOOL _nibLoaded;
    
    NSMutableDictionary *protocolSettings;
    AUIOnlineAnalystManager *analystManager;
}

/**
Implemented as [[[self alloc] init] autorelease];
 */
+ (id<AUIStimulusProtocol>)protocol;

/**
Implemented to return YES.
 */
+ (BOOL)canBeDefault;

/**
Implemented to return YES.
 */
+ (BOOL)canBeUser;

/**
Stub implementation (will throw an exception if called).
 */
- (EpochDescription*)fillNextEpoch:(EpochDescription*)epoch;

/**
Returns the protocol's Class name as protocolName.
 */
- (NSString*)protocolName;

/**
Stub implementation (will throw an exception if called).
 */
- (IBAction)showSettingsWindow:(id)sender;
/**
Returns the Class name (with " (Default)" added if this instance is a default
                        protocol).
 */
- (NSString*)settingsWindowTitle; //class name (+ (Default))


- (BOOL)validateSettings;

/**
Stub implementation (will throw an exception if called).
*/
- (void)restartUserEpochs;

/**
Stub implementation (will throw an exception if called).
 */
- (BOOL)hasNextEpoch;

/**
Enable settings in settings panel?
 Implemented as [self isDefaultProtocol] || ![self hasNextEpoch]
 */
- (BOOL)enableSettings;

/**
Several accessors are implemented (with associated instance variables).
 */
- (void)setIsDefaultProtocol:(BOOL)d;
- (BOOL)isDefaultProtocol;

- (NSMutableDictionary *)protocolSettings;
- (void)setProtocolSettings:(NSMutableDictionary *)newProtocolSettings;

- (void)setProtocolSetting:(id)newProtocolSetting forKey:(id)newKey;
- (void)removeProtocolSettingForKey:(id)newKey;

/**
Returns the Class' bundle's ID. This is probably good enough for a protocolID.
 */
+ (NSString*)protocolID;


/**
Convenience method to add a ResponseDescription object to the given epoch.
 
 @param epoch Epoch that will record the response.
 @param stream AUIDAQStreamInfo object that describes the DAQ stream to record the response from. Must be an input channel.
 @param len Duration (in seconds) to record response. Only responses from the start of the epoch can be recorded.
 @returns The response description that was added to the EpochDescription.
 */
- (ResponseDescription*)addResponseToEpoch:(EpochDescription*)epoch stream:(AUIDAQStream*)stream duration:(NSTimeInterval)len;

//stim must have parameters set.
/*!
    @method     
    @abstract   Add a stimulus to the given epoch.
    @discussion calls addStim:toEpoch:stream:setData:plot: with plot=YES.
*/

- (StimulusDescription*)addStim:(id<AUIStimulus>)stim 
                        toEpoch:(EpochDescription*)epoch 
                         stream:(AUIDAQStream*)streamInfo 
                        setData:(BOOL)setData;

/**
 @discussion Convenience method to add a StimulusDescription object to the given epoch. The stimulus'
 sampleRate will be set from the supplied AUIDAQStreamInfo, but must have its
 length and all other parameters set already.
 
 @param stim An object that conforms to the AUIStimulus protocol. All parameters + length must be set already.
 @param epoch EpochDescription that will present this stimulus.
 @param stream AUIDAQStreamInfo that describes the DAQ stream to send the stimulus out on. Must be an output channel.
 @param setData If YES, this method will generate the data from stim and set the data property of the StimulusDescription. (You probably always want to set this to YES.
 @param plot If NO, stimulus has "noplot" keyword.
 result the StimulusDescription added to epoch. If setData was NO, you must set the data yourself.
 */
- (StimulusDescription*)addStim:(id<AUIStimulus>)stim 
                        toEpoch:(EpochDescription*)epoch 
                         stream:(AUIDAQStream*)streamInfo 
                        setData:(BOOL)setData
                           plot:(BOOL)plot;

+(id)pluginType;

- (BOOL)autoExtendUserQueue;

- (NSArray*)onlineAnalysisForEpoch:(id)epoch updateDisplay:(BOOL)update;

// KVO key for settings panel. From Info.plist
- (id)settingsPanelKey;

- (void)processResponseEpoch:(Epoch*)epoch forDocument:(id)document;

- (IBAction)selectOnlineAnalysts:(id)sender;

#pragma mark Accessorizer
- (AUIOnlineAnalystManager *)analystManager;
- (void)setAnalystManager:(AUIOnlineAnalystManager *)anAnalystManager;
@end
