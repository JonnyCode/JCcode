//
//  AUIVideoStreamInfo.h
//  AcqUI
//
//  Created by Barry Wark on 8/9/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <DAQFramework/AUIIOStream.h>
#import <SpatialStimulus/SpatialStimulusPresenter.h>

@interface AUIVideoStreamInfo : AUIIOStream {
    SpatialStimulusPresenter *presenter;
    NSScreen *screen;
    
    BOOL presentationViewNibLoaded;
    IBOutlet NSOpenGLView *windowPresentationView;
    NSWindowController *presentationWindowController;
}


- (void)pushVideo:(id)renderer;

- (void)initScreen;
- (void)closeScreen;
- (void)start;
- (void)stop;
- (void)flushData;

#pragma mark Accessorizer
- (SpatialStimulusPresenter *)presenter;
- (void)setPresenter:(SpatialStimulusPresenter *)aPresenter;
- (NSScreen *)screen;
- (void)setScreen:(NSScreen *)aScreen;
@end
