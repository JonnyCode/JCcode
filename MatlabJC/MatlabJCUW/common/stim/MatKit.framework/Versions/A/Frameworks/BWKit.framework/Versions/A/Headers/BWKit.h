/*
 *  BWKit.h
 *  VizUI
 *
 *  Created by Barry Wark on 7/7/07.
 *  Copyright 2007 Physion Consulting, LLC. All rights reserved.
 *
 */

OBJC_EXPORT double BWKitVersionNumber;

#import <BWKit/BWPluginManager.h>
#import <BWKit/NSManagedObjectContext+BWAdditions.h>
#import <BWKit/NSManagedObject+BWAdditions.h>
#import <BWKit/NSFetchRequest+BWAdditions.h>
#import <BWKit/NSSortDescriptor+BWAdditions.h>
#import <BWKit/BWPluginProtocol.h>
#import <BWKit/NSDictionary+BWAdditions.h>
#import <BWKit/NSMutableDictionary+BWAdditions.h>
#import <BWKit/NSObject+BWAdditions.h>
#import <BWKit/BWPluginBase.h>
#import <BWKit/NSArray-NSSetAdditions.h>
#import <BWKit/NSMutableArray+NSSetAdditions.h>
#import <BWKit/NSString+UUID.h>
#import <BWKit/NSArray+SortAdditions.h>
#import <BWKit/NSString+TempFile.h>
#import <BWKit/NSArray+IndexingAdditions.h>
#import <BWKit/NSPredicate+BWAdditions.h>
#import <BWKit/BWAppDelegate.h>
#import <BWKit/NSException+FactoryAdditions.h>
#import <BWKit/BWErrorException.h>
#import <BWKit/NSConditionLock+BWAdditions.h>
#import <BWKit/BWBoarderlessPanel.h>
#import <BWKit/BWBoarderlessWindow.h>
#import <BWKit/BWOperation.h>
#import <BWKit/NSInvocation+BWAdditions.h>
#import <BWKit/NSMutableArray+QueueAdditions.h>
#import <BWKit/BWKitError.h>
#import <BWKit/NSError+BWAdditions.h>
#import <BWKit/BWLogger.h>
//#import <BWKit/NSObject+BWLogger.h>
#import <BWKit/BWLocalizedString.h>
#import <BWKit/BWDeprecated.h>
#import <BWKit/BWAlias.h>
#import <BWKit/NSString+RelativePathComputation.h>
#import <BWKit/BWPluginBase.h>

#import <BWKit/BWNumericData.h>
#import <BWKit/BWNumericData+TypeConversion.h>
#import <BWKit/BWNumericData+HDF5.h>
#import <BWKit/BWMutableNumericData.h>
#import <BWKit/BWNumericDataType.h>

// Core Data Entity Implementations
#import <BWKit/KeywordTag.h>
#import <BWKit/KeyValuePair.h>
#import <BWKit/NumberKeyValuePair.h>
#import <BWKit/BWFileSystemResource.h>

// NSPredicateRowTemplates for Core Data entities
#import <BWKit/BWKVPairSubqueryPredicateRowTemplate.h>

// Value Transformers
#import <BWKit/BWValueTransformerBase.h>
#import <BWKit/BWNSURLToStringTransformer.h>

// Application Kit
#import <BWKit/BWAppDelegate.h>
#import <BWKit/BWSheetController.h>
#import <BWKit/BWSingleDocumentController.h>

// QTKit
//#import <BWKit/QTMovie+CIImageFrameEnumeration.h>
//#import <BWKit/QTTrack+TimeAdditions.h>
