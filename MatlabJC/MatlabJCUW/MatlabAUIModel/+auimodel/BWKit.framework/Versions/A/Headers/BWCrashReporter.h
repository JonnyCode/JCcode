//
//  BWCrashReporter.h
//  BWKit
//
//  Created by Barry Wark on 7/10/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <Cocoa/Cocoa.h>

extern NSString * const BWCrashReporterLastCrashDateKey; //key in NSUserDefaults
extern NSString * const BWCrashReporterURLKey;

@class BWCrashReporterLog;
@class BWCrashReporter;

@protocol BWCrashReporterDelegate <NSObject>
@optional
- (BOOL)verifyResponse:(NSData*)responseData code:(NSUInteger)responseCode;
@end

@protocol BWCrashReporterSubmissionAgent <NSObject>
- (void)submitLog:(BWCrashReporterLog*)log
      userComment:(NSString*)comment 
      buildNumber:(double)buildNumber
            toURL:(NSURL*)submissionURL
   windowForSheet:(NSWindow*)theSheetWindow //non-retained; nil for application-modal
    crashReporter:(BWCrashReporter*)reporter; //non-retained
@end


/**
 Crash reporter.
 
 By default (ie when delegate is nil), uses NSUserDefaults to obtain values for submission.
 */

@interface BWCrashReporter : NSObject {
    NSString *appBundleIdentifier;
    id<BWCrashReporterDelegate> delegate;
    NSURL *submissionURL;
    NSString* version;
    
    NSWindow *submissionWindow;
    NSWindowController *windowController;
    
    NSString *userDescription;
    NSString *userSummary;
    
    id<BWCrashReporterSubmissionAgent> agent;
    
    BWCrashReporterLog *log;
}

@property (assign,readwrite) IBOutlet NSWindow *submissionWindow;
@property (retain,readwrite) NSWindowController *windowController;
@property (copy,readwrite) NSString *appBundleIdentifier;
@property (assign,readwrite) id<BWCrashReporterDelegate> delegate;
@property (retain,readonly) BWCrashReporterLog *log;

@property (copy,readwrite) NSString *userDescription;
@property (copy,readwrite) NSString *userSummary;


/**
 Designated initializer
 
 Crash logs are identified by the bundle identifier of the app that crashed.
 
 @param bundleIdentifier Bundle identifier of crashed app.
 */
- (id)initWithIdentifier:(NSString*)bundleIdentifier
           submissionURL:(NSURL*)theSubmissionURL
                 version:(NSString*)theVersion
                   agent:(id<BWCrashReporterSubmissionAgent>)theAgent;


/**
 Find and return the newest crash log that hasn't been already presented to the user.
 
 @result Latest new crash log, or nil if none exists.
 */
- (BWCrashReporterLog*)newCrashLog;

/**
 Submit a crash report.
 
 Presents a UI allowing the user to modify the report, then submits the report
 to the URL given by self.delegate or via the NSUserDefaults system.
 */
- (void)submitCrashReport:(BWCrashReporterLog*)crashLog;

- (IBAction)submit:(id)sender;
- (IBAction)cancel:(id)sender;


- (void)submissionError:(NSError*)error;
- (void)agentDidFinishSubmission:(id<BWCrashReporterSubmissionAgent>)agent;
@end