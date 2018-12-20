//
//  VUIProtocolSettingsPredicateRowTemplate.h
//  Ovation
//
//  Created by Barry Wark on 6/25/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface BWKVPairSubqueryPredicateRowTemplate : NSPredicateEditorRowTemplate <NSCoding,NSCopying> {
    NSString *collection;
    NSString *variable;
    NSDictionary *queriableProperties;
    
    NSTextField *keyView;
    NSTextField *valueView;
    NSPopUpButton *predicateOperatorTypePopup;

    NSPredicateOperatorType predicateOperatorType;
}



/*!
    @method     
    @abstract   Designated Initialzier
*/

- (id)initWithCollection:(NSString*)collectionKeyPath queriableProperties:(NSDictionary*)queriablePropertiesByName;
@end
