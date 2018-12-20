//
//  BWSheetController.h
//  BWKit
//
//  Created by Barry Wark on 2/17/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface BWSheetController : NSWindowController {
}

- (void)beginSheetModalForWindow:(NSWindow*)win
                  didEndSelector:(SEL)sel
                          target:(id)callbackTarget
                     contextInfo:(void*)contextInfo;

- (IBAction)sheetOK:(id)sender;
- (IBAction)sheetCancel:(id)sender;
@end
