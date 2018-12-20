//
//  AUIOnlineAnalystManager.h
//  AcqUI
//
//  Created by Barry Wark on 9/17/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>

extern NSString *AUIAnalystArrayChangedNotification;

@interface AUIOnlineAnalystManager : NSObject {
    IBOutlet NSPanel *analystSelectionPanel;
    IBOutlet NSArrayController *availableAnalystController;
    IBOutlet NSArrayController *analystsController;
    IBOutlet NSPanel *analystConfigPanel;
    IBOutlet NSView *configViewDestinationView;
    
    BOOL _analystSelectionNibLoaded;
    Class selectedAnalystClass; //will be a Class object
    
    NSMutableArray *onlineAnalysts;
}

@property (assign) Class selectedAnalystClass;

- (IBAction)hideAnalystSelectionPanel:(id)sender;
- (IBAction)selectOnlineAnalysts:(id)sender;
- (IBAction)addAnalyst:(id)sender;
- (IBAction)configureAnalyst:(id)sender;
- (IBAction)configureAnalystOK:(id)sender;
- (void)configAnalystSheetDidEnd:(NSWindow *)sheet returnCode:(int)returnCode contextInfo:(void *)contextInfo;
- (NSArray*)onlineAnalystClassSortDescriptors;

- (void)postArrayChangedNotification;

#pragma mark Accessorizer
- (NSMutableArray *)onlineAnalysts;
- (void)setOnlineAnalysts:(NSMutableArray *)anOnlineAnalysts;
- (void)addOnlineAnalyst:(id)anOnlineAnalyst;
- (void)removeOnlineAnalyst:(id)anOnlineAnalyst;
@end
