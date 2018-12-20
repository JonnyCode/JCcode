/*
 AUIStimulusProtocol protocol
 
 Copyright Barry Wark, 2005
 */

#import <Cocoa/Cocoa.h>

@class EpochDescription;

/*!
AUIStimulusProtocol defines the contract for interface between AcqUI and classes
 that describe experimental protocols.
 
 Roughly, a protocol should be able to provide "epochs" (defined below) in order
 to AcqUI. AcqUI will queue these ordered epochs and present them to the preparation,
 recording the response. If desired, a protocol can recieve the epoch back with its
 response included for doing post-processing or for deciding what to present next.
 
 AcqUI will queue the epochs returned from a protocol onto one of two queues: the
 "default" and the "user" queues. AcqUI will draw from the default queue
 when the user queue is empty or when the experimentor stops the user queue. Each
 queue has an associated protocol, selected by the experimentor. Thus an
 instance of a protocol may either be a default protocol--it sends epochs to the
 default queue--or a non-default protocol--it sends epochs to the user queue.
 
 Epochs that can be default epochs must be able to provide epochs quickly and without
 stop or the default queue will run dry and acquisition will stop.
 
 As a protocol writer, you should keep in mind that epochs will not be presented
 immediately. In fact, AcqUI will ask a protocol for as many epochs as it can in
 advance. Only when those epochs reach the top of their queue, and that queue is still
 running, will they be presented to the preparation. Even then there will be a small
 delay between dequeueing and presentation becaus AcqUI must buffer data to the DAQ
 interface, which has its own data buffer as well. In addition, the experimentor can
 reorder or delete epochs from either queue. Therefore, you should not assume
 that epochs that AcqUI asks your protocol for will be presented immediately, or at all.
 
 As with any protocol (or contract), you are free to define other methods for your
 protocol, to include additional UI elements or windows, to make network connections,
 or to do most anything. Keep in mind that if your protocol hogs too much CPU,
 you may prevent AcqUI from doing its job. The only way to know if your protocol takes too
 much CPU time is to try it. If you find that you want to do more processing, consider
 offloading processing to an other computer on the network, particularly 
 for post-processing work.
*/
@protocol AUIStimulusProtocol <NSObject>
/*!
@result A new instance of the protocol.
 */
+ (id<AUIStimulusProtocol>)protocol;

/*!
Protocols that can be (or make sense to be) the default protocol --
 that is the protocol that feeds AcqUI's default queue-- should
 override this method to return YES. Others should return False.
 
 @result YES if this protocol can be the default protocol.
 */
+ (BOOL)canBeDefault;

/*!
Protocols that can be (or make sense to be) the user protocol --
 that is the protocol that feeds AcqUI's experiment queue -- should
 override this method to return YES. Others should return False.
 
 @result YES if this protocol can be the user ("experiment queue") protocol.
 */
+ (BOOL)canBeUser;


/*!
 @method     
 @abstract   Validate current settings (e.g. before starting acquisition).
 @discussion AUIDocument calls validateConfig on at least the default protocol before starting acquisition. In addition, it calls validateConfig: on all protocols before calling fillNextEpoch:. Base class returns YES. Subclasses may override to do global validation of settings (you may, of course, also implement KVC-compliant validation methods for each setting (see Key-Value Coding Programming Guide: Key-Value Validation)). If the settings are not valid, returns NO and sets the appropriate error object.
 @result YES if current settigns are valid, NO otherwise. Error is set if appropriate.
 */

- (BOOL)validateSettings:(NSError**)error;


/*!
 Fills the given EpochDescription object with data for the next epoch.
 Protocol should add any stimulus and response objects, set the
 saveResponse flag and set the description property of the EpochDescription
 object. The protocol's total duration is the maximum of the duration
 of its stimuli. If the current DAQ requires data for all output channels,
 make sure that all stimuli are the same length and that all responses
 are no longer than the longest stimulus.
 
 AUIStimulusProtocolBase defines several helpful methods for adding
 stimulus and response objects to the EpochDescription.
 
 The protocol may also want to do any internal bookeeping, such as increment
 the number of epochs presented, etc.
 
 Finally, if you want, you may add info to the EpochDescription's protocolSettings
 dictionary (you may want to do this to make analysis or searching easier). You
 DO NOT need to set the Epoch's protocolSettings with the protocol's protocolSettings
 dictionary (AcqUI will take care of that for you).
 
 You may also defer this to an other protocol by calling that protocol's fillNextEpoch:.
 You may want to do this for, e.g. seal/leak testing.
 
 If you want an other protocol to post-process this epoch, set the epoch's protcol property
 to point to that protocol.
 
 @param epoch An EpochDescription object
 @result Filled EpochDescription
 */
- (EpochDescription*)fillNextEpoch:(EpochDescription*)epoch;

/*!
The user-readable name for this protocol.
 @result An NSString with a user-readable name for this protocol.
 */
- (NSString*)protocolName;

/*!
All protocols must have a settings window. The showSettingsWindow: method will
 be called by AcqUI to request that the protocol show its settings window. A
 second call when the window is already visible should close the settings window.
 */
- (IBAction)showSettingsWindow:(id)sender;

/*!
The correct title for the settings window. You may want to bind your
 settings window's title to this method to make the title visible.
 
 Title should indicate if this instance is a default protocol.
 
 E.g. 'My Protocol Settings' or 'My Protocol Settings (Default)'
 
 @result Settings window title.
 */
- (NSString*)settingsWindowTitle;


/*!
When the experimentor presses the "Go" button, they expect the current protocol
 to restart queueing epochs given its current settigns. AcqUI will send the
 current protocol a restartUserEpochs message in this case.
 */
- (void)restartUserEpochs;

/*!
AcqUI will ask the current protocol to fill epochs until it returns NO
 from hasNextEpoch. Obviously, a default epoch should never return NO
 from hasNextEpoch (AcqUI will not even call this method for the default epoch).
 
 @result Whether the protocol has more epochs to queue given current settings.
 */
- (BOOL)hasNextEpoch;

/*!
AcqUI will set the isDefaultProtocol property on a protocol instance if it
 will be queueing epochs to the default queue. Otherwise it will set the
 flag to false.
 */
- (void)setIsDefaultProtocol:(BOOL)d;
- (BOOL)isDefaultProtocol;


/*!
The protocol must provide its settings as an NSCoding-compliant dictionary.
 This dictionary will be stored in the User Defaults database and will be 
 stored in each epoch in the database to facilitate searching and analysis.
 
 @see fillNextEpoch:
 @result Dictionary of protocol's settings.
 */
- (NSMutableDictionary*)protocolSettings;
/*!
Set the settings to those provided by the given dictionary.
 
 @param d New settings.
*/
- (void)setProtocolSettings:(NSMutableDictionary*)d;

/*!
 protocolID should be a unique string identifier for this protocol.
 UTI is a good format. For example, the bandlimited white nosie protocol has ID
 
 edu.washington.bwark.acqui.protocol.WhiteNoiseStimulusProtocol
 
 The AUIStimulusProtocolBase will automatically provide the protocol's
 bundle's identifier if this method is not overriden.
 
 @result A UTI describing this stimulus.
 */
+ (NSString*)protocolID;

/*!
If YES, AcqUI will auto-extend the user epochs queue length as long as the protocol
 hasNextEpoch returns YES.
 
 Obviously, protocols that always return YES from hasNextEpoch should not return YES
 from autoExtendUserQueue.
 
@result YES if AcqUI should auto-extend the user queue length as long as this protocol
 has more epochs.
 */

- (BOOL)autoExtendUserQueue;
@end