//
//  BWLogger.h
//  BWKit
//
//  Created by Barry Wark on 2/13/08.
//  Copyright 2008 Physion Consulting, LLC. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#ifndef CHECK_FORMAT_NSSTRING
    #if MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_5
        #define CHECK_FORMAT_NSSTRING(a, b) __attribute__((format(__NSString__, a, b)))
    #else
        #define CHECK_FORMAT_NSSTRING(a, b)
    #endif
#endif

#define BWLogDebug(...)  \
[[BWLogger sharedLogger] logFuncDebug:__func__ msg:__VA_ARGS__]
#define BWLogInfo(...)   \
[[BWLogger sharedLogger] logFuncInfo:__func__ msg:__VA_ARGS__]
#define BWLogError(...)  \
[[BWLogger sharedLogger] logFuncError:__func__ msg:__VA_ARGS__]


@interface BWLogger : NSObject
{
    
}

+ (id)sharedLogger;

@end

// Helper functions that are used by the convenience GTMLogger*() macros that 
// enable the logging of function names.
// From GTMLogger.h
@interface BWLogger (BWLoggerMacroHelpers)
- (void)logFuncDebug:(const char *)func msg:(NSString *)fmt, ...
CHECK_FORMAT_NSSTRING(2, 3);
- (void)logFuncInfo:(const char *)func msg:(NSString *)fmt, ...
CHECK_FORMAT_NSSTRING(2, 3);
- (void)logFuncError:(const char *)func msg:(NSString *)fmt, ...
CHECK_FORMAT_NSSTRING(2, 3);
- (void)logFuncAssert:(const char *)func msg:(NSString *)fmt, ...
CHECK_FORMAT_NSSTRING(2, 3);
@end

//#import <asl.h>
//
//#ifdef __cplusplus 
//extern "C" {
//#endif
//void BWLogLevel(NSInteger level, id sender, char *sourceFile, NSUInteger lineNumber, NSString* format,...);
//#ifdef __cplusplus
//}
//#endif
//
//#define BWLog(format,...) BWLogLevel(ASL_LEVEL_NOTICE, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//#define BWLogInfo(format,...) BWLogLevel(ASL_LEVEL_INFO, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//#define BWLogDebug(format,...) BWLogLevel(ASL_LEVEL_DEBUG, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//#define BWLogNotice(format,...) BWLogLevel(ASL_LEVEL_NOTICE, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//#define BWLogWarning(format,...) BWLogLevel(ASL_LEVEL_WARNING, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//#define BWLogError(format,...) BWLogLevel(ASL_LEVEL_ERR, self, __FILE__, __LINE__, format, ##__VA_ARGS__)
//
//
//OBJC_EXPORT NSString * const BWLoggerExceptionName;
//
///*!
//    @class
//    @abstract    Objective-C wrapper for the ASL (Apple System Logger)
//    @discussion  Provides per-object or global loggers for the ASL. Each logger has a corresponding aslclient for each thread. Each logger stores its client in the threads' dictionary of dictionary of ASLClient wrappers (key'd by logger). When a logger needs to log from thread, it creates a new aslclient for that thread if none exists yet.
//*/
//
//@interface BWLogger : NSObject <NSCopying> {
//    NSInteger filterMask;
//    uint32_t aslClientOptions;
//    BOOL addStdErr;
//    BOOL connectImmediately;
//    NSString *facility;
//    NSFileHandle *logFile;
//}
//
//+ (BWLogger*)defaultLogger;
//+ (void)setDefaultLogger:(BWLogger*)newDefaultLogger;
//
//+ (void)setLogger:(BWLogger*)logger
//        forObject:(id)sender;
//
//+ (BWLogger*)loggerForObject:(id)obj;
//
//+ (void)clearLoggers;
//
//+ (BWLogger*)loggerWithLogURL:(NSURL*)url
//                         mask:(NSInteger)mask
//                     facility:(NSString*)facilityName
//                    addStdErr:(BOOL)stdErr
//           connectImmediately:(BOOL)connect;
//
//- (id)initWithLogURL:(NSURL*)url
//                mask:(NSInteger)mask
//            facility:(NSString*)facilityName
//           addStdErr:(BOOL)stdErr
//  connectImmediately:(BOOL)connect;
//
//- (void)closeLogFile;
//
//- (void)logLevel:(NSInteger)level
//          sender:(id)sender
//      sourceFile:(char*)sourceFile
//sourceLineNumber:(NSUInteger)lineNumber
//          format:(NSString*)format,...;
//
//- (void)logLevel:(NSInteger)level
//          sender:(id)sender
//      sourceFile:(char*)sourceFile
//sourceLineNumber:(NSUInteger)lineNumber
//          format:(NSString*)format
//      formatArgs:(va_list)formatArgs;
//
///*!
//    @method     
//    @abstract   <#(brief description)#>
//    @discussion Caller is responsible for freeing aslmsg
//*/
//
//- (aslmsg)aslMessageTemplateForSender:(id)sender
//                           sourceFile:(char*)file
//                           lineNumber:(NSUInteger)lineNumber;
//@end
