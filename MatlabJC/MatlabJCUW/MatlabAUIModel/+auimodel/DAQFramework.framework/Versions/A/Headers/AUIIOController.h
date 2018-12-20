//
//  AUIIOController.h
//  AcqUI
//
//  Created by CJ Bell on 7/10/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <AUIProtocols/AUIPluginBase.h>
#import <DAQFramework/AUIIOStream.h>
#import <DAQFramework/AUIIOStreamProperties.h>
#import <DAQFramework/AUIIOStreamIdentifier.h>

extern NSString *AUIIOControllerExceptionName;
extern NSString *AUIIOControllerErrorDomain;
extern int AUIIOControllerInvalidConfig;
extern NSString *AUIIOUnknownControllerID;
extern NSString *AUIIOControllerConfigChangedNotification;

/*!
 @protocol AUIIOController
 @abstract Protocol implemented by all I/O controllers
 @discussion  The AUIIOController protocol defines the minimum interface for I/O controllers in the DAQ Framework.
 Subclasses must call [self postConfigChangeNotification] on any config change.
 */

@protocol AUIIOController
- (NSString*)controllerID; //subclass

- (void)initHardware; //subclass NB call super if overriden
- (void)closeHardware; //subclass

- (void)startIO; //subclass. NB call after pre-loading each channel with data
- (void)stopIO; //subclass

- (void)flushQueues; //flush hardware and stream queues. subclass

- (AUIIOStream*)bathTempStream; //easy reference to the stream identified by the user as the "bath temperature" stream. all epochs must have a bath temperature. that's why bathTemp is special. subclass

    // configuration view
- (NSString*)configNibName; //subclass
- (NSView*)configView; //subclass

    //configuration settings (stored in database). subclasses must append their config to [super configDict]
    // and must call [super setConfigDict:cDict].
- (NSDictionary*)configDict;
- (void)setConfigDict:(NSDictionary*)cDict;
- (NSData*)configData; //NSCoding archived version of configDict
- (void)setConfigData:(NSData*)cData;

/*!
@method validateConfig:
 @discussion Validate the current configuration. A configuration is valid if the controller can start I/O (by calling startIO) without error.
 
 @param error NSError**. If result==NO, error will point to an allocated and autoreleased NSError describing the configuration error.
 @result YES if current config is valid. NO if not.
 */
- (BOOL)validateConfig:(NSError**)error;

- (NSSet*)streamProperties;

- (NSString *)hardwareName;
- (void)setHardwareName:(NSString *)newHardwareName;

- (NSString *)hardwareVersion;
- (void)setHardwareVersion:(NSString *)newHardwareVersion;

- (NSString *)driverVersion;
- (void)setDriverVersion:(NSString *)newDriverVersion;

- (BOOL)isRunning; //subclass
- (void)setRunning:(BOOL)flag; //subclass

- (BOOL)isOverflowed; //subclass
- (BOOL)isUnderun; //subclass

- (NSMutableSet *)channels; //this is the base storage for all streams
- (void)setChannels:(NSMutableSet *)newChannels;
- (NSArray*)channelsSortDescriptors;
- (NSArray*)channelsArray; // channels, sorted by type, then channelNumber
- (NSArray*)channelsForType:(StreamType)type includeDisabled:(BOOL)includeDisabled; //all StreamInfo's of given type (in order of channels array)
- (NSArray*)channelsForType:(StreamType)type includeDisabled:(BOOL)includeDisabled includeOwned:(BOOL)includeOwned;

- (NSArray*)analogInputChannels;
- (NSArray*)analogOutputChannels;
- (NSArray*)digitalOutputChannels;
- (NSArray*)digitalInputChannels;
- (NSArray*)videoInputChannels;
- (NSArray*)videoOutputChannels;
- (NSArray*)outputChannels;
- (NSArray*)inputChannels;
- (NSArray*)enabledChannels;

    // these retun only non-owned (by external device) channels
- (NSArray*)availableAnalogInputChannels;
- (NSArray*)availableAnalogOutputChannels;
- (NSArray*)availableDigitalOutputChannels;
- (NSArray*)availableDigitalInputChannels;
- (NSArray*)availableVideoInputChannels;
- (NSArray*)availableVideoOutputChannels;

- (BOOL)hardwareInitialized; //subclass
- (void)setHardwareInitialized:(BOOL)flag; //subclass

    // returns nil if no stream found
- (AUIIOStream*)streamInfoForChannel:(unsigned short)channel 
                              ofType:(StreamType)type; //depreciated
- (AUIIOStream*)streamForChannel:(unsigned short)channel
                          ofType:(StreamType)type;

    /*
     returns a dictionary of the properties for given stream (as read from config data).
     Stream properties are defined by the stream's streamPropertyKeys. Since stream properties are handled by AUIIOController directly, it can be done as a class method (no specific controller instance required).
     */
+ (NSDictionary*)streamPropertiesForChannel:(unsigned short)channelNumber
                                     ofType:(StreamType)type
                                 configData:(NSData*)configData;

+ (NSDictionary*)streamPropertiesForStreamIdentifier:(AUIIOStreamIdentifier*)streamIdentifier
                                          configData:(NSData*)configData;

- (BOOL)canParseConfigDataForID:(NSString*)daqID;

    //subclass
-(void)writeToFIFO:(unsigned short)channel 
            ofType:(StreamType)type 
              data:(id)data;

    //subclass. throws AUIIOException if aquisition fails or timing is wrong (takes to short/long)
-(id)readFromFIFO:(unsigned short)channel 
           ofType:(StreamType)type 
       numSamples:(unsigned long)nsamples 
      readStreams:(NSSet*)rStreams;
@end

@interface AUIIOController : AUIPluginBase <AUIIOController> {
    NSString *hardwareName;
	NSString *hardwareVersion;
	NSString *driverVersion;
	NSMutableSet *channels;
    
    BOOL postConfigChangeNotificationsForStreamChanges;
    BOOL hasSuppressedConfigChangeNotifications;
}

@property (assign,readwrite) BOOL postConfigChangeNotificationsForStreamChanges;

+(NSArray*)registeredControllers;
+(void)registerController:(AUIIOController*)c;
+(AUIIOController*)registeredControllerForID:(NSString*)cID;
+(AUIIOController*)registeredControllerToParseConfigDataForID:(NSString*)cID;
+(NSArray*)availableControllers; //interfaces for this particular controller type

+(void)registerControllerClass:(Class) newControllerClass; //call to register a class to be initialized (e.g. call from +initialize)
+(void)initializeControllers:(id)delegate; //to be called by e.g. AppDelegate to initialize all registered controller classes. delegate may be nil.
+(void)setup; //will be called by initializeControllers for each registered controller class. subclass

/*!
@method     currentController
 @abstract   Global current I/O controller
 @discussion AUIIOController maintains a global current controller. Returns the current controller (or nil) if none is set. The client application can make use of this parameter as it wants. The current controller is used by the DAQ Framework only to help in unarchiving AUIIOStreamInfos. If the controller that owned a stream is not registered when the stream is unarchived, and the current controller can unarchive data from that (missing) controller, the corresponding stream of the current controller is substituted for the arhived stream.
 @seealso AUIIOStreamReference
 */
+ (id<AUIIOController>)currentController;
/*!
    @method     setCurrentController
    @abstract   Set the global current controller for the DAQ Framework.
    @discussion (comprehensive description)
*/

+ (void)setCurrentController:(id<AUIIOController>)newCurrentController;

// NB: Any method may throw AUIIOException (or a subclass)

/*
 controllerID may not change after initialization, since it's used to hash AUIIOStreamInfo's that have this controller
 as their controller
 */
- (NSString*)controllerID; //subclass

- (id)initWithNumChannels:(unsigned)nChannels initHardware:(BOOL)initHardware; //subclass
- (void)initHardware; //subclass NB call super if overriden
- (void)closeHardware; //subclass

- (void)startIO; //subclass. NB call after pre-loading each channel with data
- (void)stopIO; //subclass

- (void)flushQueues; //flush hardware and stream queues. subclass

- (AUIIOStream*)bathTempStream; //easy reference to the stream identified by the user as the "bath temperature" stream. all epochs must have a bath temperature. that's why bathTemp is special. subclass

// configuration view
- (NSString*)configNibName; //subclass
- (NSView*)configView; //subclass

//configuration settings (stored in database). subclasses must append their config to [super configDict]
// and must call [super setConfigDict:cDict].
- (NSDictionary*)configDict;
- (void)setConfigDict:(NSDictionary*)cDict;
- (NSData*)configData; //NSCoding archived version of configDict
- (void)setConfigData:(NSData*)cData;

/*!
 @method
 @discussion Verifies that all analog streams have a valid external device. Subclasses may override this method (but should call super's implementation) to do further validation of their configuration. Called by AUIDocument before starting queues.
 */
- (BOOL)validateConfig:(NSError**)error;

- (NSSet*)streamProperties;

- (NSString *)hardwareName;
- (void)setHardwareName:(NSString *)newHardwareName;

- (NSString *)hardwareVersion;
- (void)setHardwareVersion:(NSString *)newHardwareVersion;

- (NSString *)driverVersion;
- (void)setDriverVersion:(NSString *)newDriverVersion;

- (BOOL)isRunning; //subclass
- (void)setRunning:(BOOL)flag; //subclass

- (BOOL)isOverflowed; //subclass
- (BOOL)isUnderun; //subclass

- (NSMutableSet *)channels; //this is the base storage for all streams
- (void)setChannels:(NSMutableSet *)newChannels;
- (NSArray*)channelsSortDescriptors;
- (NSArray*)channelsArray; // channels, sorted by type, then channelNumber
- (NSArray*)channelsForType:(StreamType)type includeDisabled:(BOOL)includeDisabled; //all StreamInfo's of given type (in order of channels array)
- (NSArray*)channelsForType:(StreamType)type includeDisabled:(BOOL)includeDisabled includeOwned:(BOOL)includeOwned;

- (NSArray*)analogInputChannels;
- (NSArray*)analogOutputChannels;
- (NSArray*)digitalOutputChannels;
- (NSArray*)digitalInputChannels;
- (NSArray*)videoInputChannels;
- (NSArray*)videoOutputChannels;
- (NSArray*)outputChannels;
- (NSArray*)inputChannels;
- (NSArray*)enabledChannels;

// these retun only non-owned (by external device) channels
- (NSArray*)availableAnalogInputChannels;
- (NSArray*)availableAnalogOutputChannels;
- (NSArray*)availableDigitalOutputChannels;
- (NSArray*)availableDigitalInputChannels;
- (NSArray*)availableVideoInputChannels;
- (NSArray*)availableVideoOutputChannels;

- (BOOL)hardwareInitialized; //subclass
- (void)setHardwareInitialized:(BOOL)flag; //subclass

// returns nil if no stream found
- (AUIIOStream*)streamInfoForChannel:(unsigned short)channel 
                              ofType:(StreamType)type;
- (AUIIOStream*)streamForChannel:(unsigned short)channel
                          ofType:(StreamType)type;

/*!
 @abstract Get the dictionary of stream properties for a specified stream in a DAQ config data.
 @discussion Stream properties are defined by the stream's streamPropertyKeys
 @result A dictionary of the properties for given stream (as read from config data) or nil if there is no AUIIOStreamProperties in config data with the given type and channelNumber.

 */
+ (NSDictionary*)streamPropertiesForChannel:(unsigned short)channelNumber
                                     ofType:(StreamType)type
                                 configData:(NSData*)configData;

+ (NSDictionary*)streamPropertiesForStreamIdentifier:(AUIIOStreamIdentifier*)streamIdentifier
                                          configData:(NSData*)configData;

/*!
    @method     
    @abstract   <#(brief description)#>
    @discussion <#(comprehensive description)#>
 @throws NSException if streamIdentifier.controllerID does not match self.controllerID
*/

- (AUIIOStream*)streamForStreamIdentifier:(AUIIOStreamIdentifier*)streamIdentifier;

- (NSPredicate*)channelIsEnabledPredicate;

-(void)appWillQuit:(NSNotification*)n;
- (void)handleStreamChangedNotification:(NSNotification*)n;

- (BOOL)canParseConfigDataForID:(NSString*)daqID;

//subclass
-(void)writeToFIFO:(unsigned short)channel 
            ofType:(StreamType)type 
              data:(id)data;

//subclass. throws AUIIOException if aquisition fails or timing is wrong (takes to short/long)
-(id)readFromFIFO:(unsigned short)channel 
           ofType:(StreamType)type 
       numSamples:(unsigned long)nsamples 
      readStreams:(NSSet*)rStreams;


/*!
    @method     
    @abstract   Allows sender to turn of configChangeNotifications in response to streamChangeNotifications from self.streams
    @discussion If set to NO, will post a single configChangeNotification when later set to YES, if any would have been posted in the mean time.
*/

- (void)setPostConfigChangeNotificationsForStreamChanges:(BOOL)post;

/*!
    @method     
    @abstract   Private method
    @discussion Subclasses should call when config changes to post ConfigChangedNotification.
*/

- (void)postConfigChangedNotification;
@end
