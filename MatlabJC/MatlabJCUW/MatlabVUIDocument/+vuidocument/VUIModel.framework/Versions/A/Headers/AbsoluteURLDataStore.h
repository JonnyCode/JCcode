//
//  AbsoluteURLDataStore.h
//  Ovation
//
//  Created by Barry Wark on 2/16/08.
//  Copyright 2008 Barry Wark. All rights reserved.
//

#import <CoreData/CoreData.h>
#import "URLDataStore.h"


@interface AbsoluteURLDataStore :  URLDataStore  
{
}

-(void)setUrl:(NSURL*)theURL;
@end


