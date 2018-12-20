/*
 *  AUIModel.h
 *  AcqUI
 *
 *  Created by Barry Wark on 12/8/07.
 *  Copyright 2007 Barry Wark. All rights reserved.
 *
 */

extern NSString * const AUIModelErrorDomain;
extern const NSInteger AUIModelExceptionError;
extern const NSInteger AUIModelDatabaseLogicError;
extern const NSInteger AUIModelDataMigrationError;

#import <AUIModel/AUIModelFunctions.h>

#import <AUIModel/Experiment.h>
#import <AUIModel/RecordedCell.h>
#import <AUIModel/Epoch.h>
#import <AUIModel/IOBase.h>
#import <AUIModel/Stimulus.h>
#import <AUIModel/Response.h>
#import <AUIModel/Note.h>
#import <BWKit/NumberKeyValuePair.h>
#import <AUIModel/AUIKeyValuePair.h>
#import <AUIModel/DAQConfigContainer.h>
#import <AUIModel/StreamKeyValuePair.h>
#import <BWKit/KeywordTag.h>

#import <AUIModel/RecordedCellAdditions.h>
#import <AUIModel/StimulusAdditions.h>
#import <AUIModel/Epoch-Additions.h>
#import <AUIModel/AUIKeyValuePair.h>
#import <AUIModel/IOBaseDAQAdditions.h>
#import <AUIModel/ResponseAdditions.h>
#import <AUIModel/StreamKeyValuePairAdditions.h>
#import <AUIModel/AUIFileSystemResource.h>

#import <AUIModel/AUIModelMigrator.h>

#import <AUIModel/CachedData.h>