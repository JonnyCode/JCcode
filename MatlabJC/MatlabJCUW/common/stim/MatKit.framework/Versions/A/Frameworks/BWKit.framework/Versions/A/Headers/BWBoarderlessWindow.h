//
//  BWBoarderlessWindow.h
//  AcqUI
//
//  Created by Barry Wark on 8/13/06.
//  Copyright 2006 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface BWBoarderlessWindow : NSWindow {

}

- (id)initWithContentRect:(NSRect)contentRect styleMask:(NSUInteger)styleMask backing:(NSBackingStoreType)bufferingType defer:(BOOL)deferCreation;

- (BOOL)canBecomeKeyWindow;
@end
