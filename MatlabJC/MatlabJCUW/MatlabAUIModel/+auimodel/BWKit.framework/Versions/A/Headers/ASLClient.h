//
//  ASLClient.h
//  BWKit
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <asl.h>

@interface ASLClient : NSObject {
    aslclient client;
}

@property (assign,readwrite) aslclient client;

+ (ASLClient*)wrapperForASLClient:(aslclient)aslclient;
- (id)initWithClient:(aslclient)aslClient;
@end
