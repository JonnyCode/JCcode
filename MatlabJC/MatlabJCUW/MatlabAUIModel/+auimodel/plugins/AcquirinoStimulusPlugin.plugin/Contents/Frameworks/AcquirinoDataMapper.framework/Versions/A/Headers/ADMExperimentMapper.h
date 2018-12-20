//
//  ADMExperimentMapper.h
//  AcquirinoDataMapper
//
//  Created by Daniel Carleton on 1/7/08.
//  Copyright 2008 __MyCompanyName__. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <CoreData/CoreData.h>

#import "ADMStimulusParameterMapper.h"
#import "ADMDataFileReader.h"

extern NSString* const ADMAcquirinoStimulusVersionKey;
extern NSString* const ADMAcquirinoStimulusSamplingIntervalKey;
extern NSString* const ADMAcquirinoStimulusDataPointsKey;

@interface ADMExperimentMapper : NSObject {
@private
  NSString* notebookFilePath;
  NSSet* dataFilePaths;
  NSData* daqConfig;
  NSManagedObjectContext* context;
  NSURL* auisqlUrl;
  NSDate* experimentStartDate;
  ADMStimulusParameterMapper* stimulusMapper;
  NSSet* excludeFromRollupSettings;
  NSURL* hdf5FileUrl;
  ADMItcDevice itcDevice;
}

@property (readonly,copy) NSString* notebookFilePath;
@property (readonly,assign) NSSet* dataFilePaths;
@property (readonly,assign) NSData* daqConfig;
@property (readonly,assign) NSManagedObjectContext* context;
@property (assign) NSURL* auisqlUrl;
@property (assign) NSURL* hdf5FileUrl;

+ (id) mapperForNotebook: (NSString *) notebookPath
            andDataFiles: (NSSet *) cellPaths
           withDaqConfig: (NSData *) daqConfigData
             intoContext: (NSManagedObjectContext *) moc
         havingAuisqlUrl: (NSURL *) auisqlUrl
             usingDevice: (ADMItcDevice) device;

- (id) initWithNotebook: (NSString *) notebookPath
           andDataFiles: (NSSet *) cellPaths
          withDaqConfig: (NSData *) daqConfigData
            intoContext: (NSManagedObjectContext *) moc
        havingAuisqlUrl: (NSURL* ) auisqlUrl
            usingDevice: (ADMItcDevice) device;

- (void) map;

@end

