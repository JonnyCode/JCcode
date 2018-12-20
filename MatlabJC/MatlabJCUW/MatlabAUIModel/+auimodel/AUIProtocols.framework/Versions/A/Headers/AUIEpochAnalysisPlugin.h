/*
 *  AUIOnlineAnalysisPlugin.h
 *  AcqUI
 *
 *  Created by Barry Wark on 6/11/07.
 *  Copyright 2007 Barry Wark. All rights reserved.
 *
 */

#import <AUIProtocols/AUIPluginProtocol.h>

#define AUI_EPOCH_ANALYSIS_PLUGIN_TYPE @"edu.washington.bwark.acqui.epoch_analyst"

@protocol AUIEpochAnalysisPlugin <AUIPlugin,NSCoding>
//create and return an initialized, autoreleased analyst
+ (id<AUIEpochAnalysisPlugin>)analyst;
+ (BOOL)userVisible; //some may want to be protocol-callable only (e.g. R&C) via the protocol property of EpochDescription

+ (NSString*)userName; //name of analyst type
- (NSString*)uiShortDescription; //title of this analyst

//analyze given epoch (type AUIModel/Epoch)
- (id)onlineAnalysisForEpoch:(id)epoch updateDisplay:(BOOL)update;
- (id)onlineAnalysisForEpochs:(NSArray*)epochs updateDisplay:(BOOL)update;

- (BOOL)enabled;
- (void)setEnabled:(BOOL)newEnabled;

- (NSView*)configView; //app will handle display of this config view (probably in a sheet). should autosize.
- (NSView*)displayView; //app will handle display of this view, probably in a utility window. should autosize.

- (IBAction)reset:(id)sender; //reset analysis

@optional
- (id)offlineAnalysisForIOBases:(NSArray*)ioBases updateDisplay:(BOOL)update; //run analysis, updating displayView if requested. ioBases are provided to allow analysis on diverse streams etc.
@end