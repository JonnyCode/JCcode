//
//  URLDataStore.h
//  Ovation
//
//  Created by Barry Wark on 2/17/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import "DataStore.h"


@interface URLDataStore :  DataStore  
{
}

@property (retain) NSURL* url;
@end


