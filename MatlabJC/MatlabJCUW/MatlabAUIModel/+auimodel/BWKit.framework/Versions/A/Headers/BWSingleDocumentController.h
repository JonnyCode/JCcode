//
//  BWSingleDocumentController.h
//  BWKit
//
//  Created by Barry Wark on 5/19/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>


/*!
    @class
    @abstract    Overrides NSDocumentController behavior to only permit 
 one open document.
    @discussion  The application delegate should also override 
 -application:openFile: to call the document controller's 
 -openDocumentWithContentsOfFile:display: if there are no open documents
 (returning YES in either case).
*/

@interface BWSingleDocumentController : NSDocumentController {

}

- (void)openDocument:(id)sender;
- (void)newDocument:(id)sender;
- (NSDocument*)document; //the single document, or nil if none
- (BOOL)validateUserInterfaceItem:(id < NSValidatedUserInterfaceItem >)anItem;
@end
