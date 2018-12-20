//
//  SpotlightDataStore.h
//  Ovation
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import "DataStore.h"


@interface SpotlightDataStore :  DataStore  
{
}

@property (retain) NSString * queryString;

@end


