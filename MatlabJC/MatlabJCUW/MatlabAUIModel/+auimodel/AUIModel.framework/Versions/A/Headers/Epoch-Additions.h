//
//  Epoch-Additions.h
//  AcqUI
//
//  Created by Barry Wark on 1/30/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <AUIModel/Epoch.h>
#import <DAQFramework/AUIDAQStream.h>
#import <AUIModel/AUIKeyValuePair.h>
#import <BWKit/BWErrorException.h>

@class Stimulus;
@class EpochDescription;

@interface Epoch (EpochAdditions)


- (id)stimulusWithChannelID:(int64_t)channelID ofType:(StreamType)type;
- (id)responseWithChannelID:(int64_t)channelID ofType:(StreamType)type;
- (id)ioBaseWithChannelID:(int64_t)channelID ofType:(StreamType)type;
- (NSArray*)ioBasesWithChannelID:(int64_t)channelID ofType:(StreamType)type;

/* 
    These return an NSArray of the given IOBase type with a streamInfo whose user description matches userDescription
*/
- (NSArray*)stimuliWithStreamName:(NSString*)userDescription;
- (NSArray*)responsesWithStreamName:(NSString*)userDescription;

/*!
 @method     
 @abstract   Find IOBases from epochs.
 @discussion Finds all IOBases of given type/userDescription from epochs. Returned array is sorted
 by sortKey (a key on Epochs).
 @param streamUserDescription The userDescription of the IOBases' AUIIOStream.
 @param sortKey Key on Epoch to sort the resulting array. If nil, ioBases are in order of `epochs`.
*/

+ (NSArray*)ioBasesForEpochs:(NSArray*)epochs 
                    metaType:(StreamMetaType)type 
             userDescription:(NSString*)streamUserDescription
                     sortKey:(NSString*)sortKey;

/*!
 @method     
 @abstract   Find IOBases from Epochs, sorted by Epoch.startDate
 @discussion Convenience method. Calls ioBasesForEpochs:metaType:userDescription:sortKey: with sortKey nil.
*/

+ (NSArray*)ioBasesForEpochs:(NSArray*)epochs 
                    metaType:(StreamMetaType)type 
             userDescription:(NSString*)streamUserDescription;

/*!
 @method     
 @abstract   Find IOBases from epochs.
 @discussion Finds all IOBases of given channelID/type from epochs. Returned array is sorted
 by sortKey (a key on Epochs).
 @param streamUserDescription The userDescription of the IOBases' AUIIOStream.
 @param sortKey Key on Epoch to sort the resulting array. If nil, ioBases are in order of `epochs`.
 */
+ (NSArray*)ioBasesForEpochs:(NSArray*)epochs 
                   channelID:(int64_t)channelID 
                      ofType:(StreamType)type
                     sortKey:(NSString*)sortKey;

/*!
 @method     
 @abstract   Find IOBases from Epochs, sorted by Epoch.startDate
 @discussion Convenience method. Calls ioBasesForEpochs:channelID:type:sortKey: with sortKey nil.
 */
+ (NSArray*)ioBasesForEpochs:(NSArray*)epochs 
                   channelID:(int64_t)channelID 
                      ofType:(StreamType)type;

- (NSDictionary*)protocolSettings;
- (void)setProtocolSettings:(NSDictionary*)d;
@end
