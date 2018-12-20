//
//  KeyValuePair.h
//  AcqUI
//
//  Created by Barry Wark on 11/7/07.
//  Copyright 2007 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>


OBJC_EXPORT NSString * const KVPValueTypeKey;
OBJC_EXPORT NSString * const KVPValueRequiredProtocolKey;
OBJC_EXPORT NSString * const KVPValueKeyKey;
OBJC_EXPORT NSString * const KVPManagedObjectClassKey;

/**
    @protocol
    @brief    Formal protocol defining methods required of KeyValuePair subclasses
*/

@protocol KeyValuePairSubEntityMethods
/**
 @brief   Return a (sub)entity to store value
 
 Sub-entities of KeyValuePair will have this method called when KeyValuePair determines that they are the correct entity for storing value. Subclasses may use this method to further refine the entity used to store value (i.e. to implement a class cluster).
 */

+ (NSEntityDescription*)entitySubclassForValue:(NSValue*)value
                        inManagedObjectContext:(NSManagedObjectContext*)moc;
@end


/**
 @class
 @brief    	Abstract entity for storing key-value pairs in a CoreData store.
 
 KeyValuePair is an abstract entity, analogous to a class cluster in Cocoa (e.g. NSNumber). Sub-entities store value sof specific types. You should not work directly with KeyValuePair's sub-entities. Instead, use the factory methods defined below to create a KeyValuePair sub-entity. All sub-entities have a value property which holds the value of that sub-entity's type.
 
 Sub-entities define a list of object classes that they can store, along with (optionally) a protocol that stored objects must conform to (e.g. NSCoding). The list of storeable classes and required protocol(s) is stored in the sub-entitie's userInfo dictionary in the managed object model. List of classes should be from most specific to least (e.g. AUIDAQStream AUIIOStream NSObject) and be white-space delimited. Only entities advertising the most specific (e.g. closest in ancestry to [value class]) class are considered. If more than one entity can store a value, an NSGenericException is raised with format stating the entities (and value class). Thus there can be only one KeyValuePair subentity that advertises the class/protocol combination closest to [value class] of all the class/protocol advertised.
 
 Sub-entities may further refine the (sub)entity used for values that they store by implementing the KeyValuePairSubclassMethods informal protocol (entitySubclassForValue:inManagedObjectContext:). These entities should additionally define KVPManageObjectClass => [entities managed object class name] in their entity's userInfo dictionary in case access to the entitySubclassForValue:inManagedObjectContext: method is needed during data migration.
 
 All incoming relationships on a KeyValuePair (i.e. those representing sets of KVPs) should be set to "Deny" delete rule to prevent deletion of a KVP during set update if other sets contain that KVP. Return connections should be marked "to-many".
 @updated 1-8-2009
 */

@interface KeyValuePair :  NSManagedObject  
{
}

@property (retain) NSNumber * boolValue;
@property (retain) id codedValue;
@property (retain) NSData * dataValue;
@property (retain) NSDate * dateValue;
@property (retain) NSNumber * doubleValue;
@property (retain) NSNumber * floatValue;
@property (retain) NSNumber * intValue;
@property (retain) NSString * key;
@property (retain) NSNumber * longIntValue;
@property (retain) NSString * stringValue;
@property (readonly) id value;


/**
 @method
 @brief    Must be overriden to return the managed object model that this KeyValuePair entity is defined in.
 
 @throw NSNotImplementedException
*/

+ (NSManagedObjectModel*)managedObjectModel;


/**
 @method     
 @brief   All property names that may contain a value (i.e. not 'key')
 
 All properties with name xxxValue.
 */

+ (NSSet*)valuePropertyNames;

/**
 @method     
 @brief   All attribute value property names that can be queried in a predicate (no data values; no relationship values)
 
 Excludes data values.
 */

+ (NSSet*)queriableValuePropertyNamesInManagedObjectModel:(NSManagedObjectModel*)model;
+ (NSDictionary*)queriableValuePropertiesByNameInManagedObjectModel:(NSManagedObjectModel*)model;


/**
 @brief   Find unique KVP in set
 
 Checks that KVP with given key is unique in set, and returns it, if found.
 @param key KVP key to search for.
 @param kvpSet Set of KVPs to search within.
 @throw NSGenericException if key is not unique in kvpSet.
 @updated 12-17-2007
 */

+ (KeyValuePair*)checkedUniqueKVPWithKey:(NSString*)key
                                   inSet:(NSSet*)kvpSet;
/**
 @brief   	Factory method for insertion of KeyValuePair sub-entities.
 
 This factory method is the designated way to get or create a KVP. Like insertNewObjectForEntityForName:inManagedObjectContext:, this method returns an autoreleased NSManagedObject. Analogous to factory methods for other Cocoa class clusters, this method returns an NSManagedObject whose entity is a sub-entity of KeyValuePair. The particular sub-entity is chosen by the list of value types each sub-entity declares that it can store. For example, an NSString value will be stored in a StringKeyValuePair. See KeyValuePair(Additions) documentation for more info.
 
 @throw		NSGenericException if [value class] is not a storeable type.
 @throw		NSGenericException if more than one entity advertises [value class] as one of its storable classes.
 @upadted 12-17-2007
 */
+ (id)insertNewKeyValuePairWithKey:(NSString*)key 
                             value:(id)value 
            inManagedObjectContext:(NSManagedObjectContext*)moc;


/**
 @brief   Build a dictionary of key-value pairs from a set of KeyValuePairs.
 
 Returned dictionary is NOT a proxy for the original keys/values. Changes made to result will not directly affect the KVP set. Use updateKVPSet:owner:fromDictionary: if you want to make changes.
 @param kvps Set of KeyValuePair subentities.
 @return Dictionary of key-value pairs corresponding to the members of kvps.
 @updated 12-17-2007
 */

+ (NSDictionary*)dictionaryFromKeyValuePairs:(NSSet*)kvps;

/**
 @brief   Updates a mutable set of KeyValuePairs from a dictionary of key-value pairs
 
 KeyValuePairs are created in owners.managedObjectContext. KeyValuePairs whose key does not appear ni d will be deleted from kvpSet (and from their managedObjectContext). KeyValuePairs whose key appears in d but whose value is different from the corresponding value in d will either have their values updated (if the type of the new value is still storeable in the original KVP) or deleted and recreated if a new value type is given in d.
 @param kvpSet NSMutableSet of KeyValuePairs
 @param	d Dictionary of key-value pairs. If d == nil, this method assumes empty dictionary.
 @param owner Owner of kvpSet. Defines managedObjectContext in which to create new KVPs.
 @throw NSGenericException if more than one entity advertises any d.value.class as one of its storable classes.
 @depreciated
 @updated 12-17-2007
 */

+ (void)updateKVPSet:(NSMutableSet*)kvpSet
               owner:(NSManagedObject*)owner
      fromDictionary:(NSDictionary*)d;

/**
 @brief   Updates a mutable set of KeyValuePairs from a dictionary of key-value pairs
 KeyValuePairs are created in moc. KeyValuePairs whose key does not appear ni d will be deleted from kvpSet (and from their managedObjectContext). KeyValuePairs whose key appears in d but whose value is different from the corresponding value in d will either have their values updated (if the type of the new value is still storeable in the original KVP) or deleted and recreated if a new value type is given in d.
 @param kvpSet NSMutableSet of KeyValuePairs
 @param	d Dictionary of key-value pairs. If d == nil, this method assumes empty dictionary.
 @param moc NSManagedObjectContext* to create new KVPs in.
 @throw NSGenericException if more than one entity advertises any d.value.class as one of its storable classes.
 
 @updated 12-17-2007
 */

+ (void)updateKVPSet:(NSMutableSet*)kvpSet
      fromDictionary:(NSDictionary*)d
managedObjectContext:(NSManagedObjectContext*)moc;

/**
 @brief   Access kvp value
 The property storing kvp.value is dependent of the value's type. Using this method, you can access the value without knowing its type.
 @param kvp A KeyValuePair sub-entity
 @return kvp's value.
 */

+ (id)kvpValueForInstance:(NSManagedObject*)kvp;


/**
 @brief   Access instance's value
 Sub-entities which define subclasses of KeyValuePair can use this method to access instances' value directly.
 @see kvpValueForInstance:
 */
- (id)value;



/**
 Generate a predicate suitable for quering KeyValue entities.
 
 Generates a SUBQUERY predicate from given collection for KeyValue pairs of given key and comparison value with comparison predicate modifier and type.
 
 @param collectionKeyPath Collection key path of SUBQUERY
 @param variable SUBQUERY variable name
 @param key KeyValue.key (i.e. kvp.key == key)
 @param value Comparison value (e.g. kvp.value == value)
 @param modifier Comparison predicate modifier
 @param type Predicate operator type
 @param queriableProperties Result of desired subclass' +queriablePropertiesByName
 @return NSPredicate suitable for quering KeyValue entities.
 */
+ (NSPredicate*)predicateForCollection:(NSString*)collectionKeyPath
                              variable:(NSString*)variable
                                    key:(NSString*)key
                                 value:(id)value
                              modifier:(NSComparisonPredicateModifier)modifier
                                  type:(NSPredicateOperatorType)type 
             queriablePropertiesByName:(NSDictionary*)queriableProperties;

+ (NSString*)entityNameForValue:(id)value 
           managedObjectContext:(NSManagedObjectContext*)moc;
@end
