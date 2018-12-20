
@class Epoch;

@interface NSObject (AUIStimulusProtocolAdditions)
-(void)processResponseEpoch:(Epoch*)currEpoch forDocument:(id)document;
@end