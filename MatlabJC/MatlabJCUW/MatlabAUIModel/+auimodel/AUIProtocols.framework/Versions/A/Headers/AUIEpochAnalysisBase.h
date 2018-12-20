//
//  AUIOnlineAnalysisBase.h
//  AcqUI
//
//  Created by Barry Wark on 6/11/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "AUIPluginBase.h"
#import "AUIEpochAnalysisPlugin.h"

@interface AUIEpochAnalysisBase : AUIPluginBase <AUIEpochAnalysisPlugin> {
    BOOL _enabled;
    BOOL _nibLoaded;
    IBOutlet NSView *configView;
    IBOutlet NSView *displayView;
}
+ (id<AUIEpochAnalysisPlugin>)analyst;

+ (BOOL)userVisible; //some may want to be protocol-callable only (e.g. R&C) via the protocol property of EpochDescription
+ (NSString*)userName;

- (NSString*)uiShortDescription;
    //analyze given epoch (type AUIModel/Epoch)
- (id)onlineAnalysisForEpoch:(id)epoch updateDisplay:(BOOL)update;;
- (id)onlineAnalysisForEpochs:(NSArray*)epochs updateDisplay:(BOOL)update;; //should show GUI/plots

- (BOOL)enabled;
- (void)setEnabled:(BOOL)newEnabled;

- (NSView*)configView;
- (NSView*)displayView;
+ (NSString*)nibName; //must override if it's not <ClassName>.nib

- (void)encodeWithCoder:(NSCoder*)encoder;
- (id)initWithCoder:(NSCoder*)decoder;

- (IBAction)reset:(id)sender;
@end