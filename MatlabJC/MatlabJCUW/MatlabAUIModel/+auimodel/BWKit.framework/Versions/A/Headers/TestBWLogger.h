//
//  TestBWLogger.h
//  BWKit
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <SenTestingKit/SenTestingKit.h>
#import <asl.h>


@interface TestBWLogger : SenTestCase <NSCopying> {
    id client1;
    id client2;
    NSConditionLock *lock;
}

- (void)testBWLoggerLogsFromSeparateClientForEachThread;
- (void)testBWLoggerChoosesCorrectLoggerForSender;
- (void)testBWLoggerLogsToASL;
- (void)testBWLoggerFiltersLevels;
- (void)testBWLoggerLogsToURL;

- (void)logFromSeparateThread;
@end
