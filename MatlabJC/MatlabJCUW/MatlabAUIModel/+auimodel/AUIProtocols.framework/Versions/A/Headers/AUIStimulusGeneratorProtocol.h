/**
 AUIStimulus protocol
 
 AUIStimulus defines the methods of a stimulus generator for the AcqUI
 application.
 
 Copyright Barry Wark, 2005
 */

@protocol AUIStimulus <NSObject>

@property (assign,readwrite) NSTimeInterval length;
@property (assign,readwrite) double sampleRate;

/**
Logically equivalent to 
 id<AUIStimulus> stim;
 
 [stim stimulusDataForParamsDictionary:[self paramsDictionary]];
 
 @returns A data object with the stimulus for the current parameters.
 */
-(NSData*)stimulusData;

/**
The returned data should be of self.length at self.sampleRate and should be
 in the natural units. AcqUI will take care of converting to e.g., DAQ counts etc.
 
 Calls stimulusDataForParamsDictionary:p version:[self version] fileSystemResources:nil.
 */
-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p;

/**
Same as -stimulusDataForParamsDictionary:, but computes stimulus as given version.
 
 PRECONDITION: version<=[self version]
 */
-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p version:(NSUInteger)version;

/**
 Designated method to reconstruct data.
 
 @param p A dictionary of parameters (specific keys are defined by the stimulus class)
 @param version Stimulus generation code version (defined by stimulus plugin).
 @param fileSystemResources Set of file system resource URLs associated with a stimulus. May be nil.
 @param subStimuli Set of BWNumericData* for each sub stimulus. May be nil.
 @return A data object with the stimulus for the parameters contained in the argument.
 */
-(NSData*)stimulusDataForParamsDictionary:(NSDictionary*)p version:(NSUInteger)version fileSystemResources:(NSSet*)fileSystemResources subStimuli:(NSSet*)subStimuli;

/**
 Dictionary defining stimulus parameters. A stimulus must be able to reconstruct
 a stimulus given this dictionary.
 
 Since this dictionary might be re-used by AcqUI, the stimulus generator must
 not change its state (i.e. re-seed is RNG) unless told to by an other object.
 
 Stimuli that use outside libraries to generate their stimulus should
 include those libraries' versions in their parameters.
 
 You do not need to include the stimulus ID, length or sampleRate.
 
 @returns Dictionary of stimulus parameters.
 */
-(NSDictionary*)paramsDictionary;


/*
 stimulusID should be a unique string identifier for this stimulus type.
 UTI is a good format. For example, the step stimulus is defined as
 
 edu.washington.bwark.acqui.stimulus.step
 
 The stimulus ID does not need to include the stimulus version, but should
 be unique for each stimulus type (probably equivalent to classes).
 
 @returns A UTI describing this stimulus.
 */
 
+(NSString*)stimulusID;

/*
 The version of this implementation. If possible, this version will be used
 to reconstruct the stimulus in the future, e.g. for analysis. If you make
 changes to the code for stimulus generation that would affect how previously
 generated stimuli will be reconstructed, or if you change the stimulus
 parameters, you should increment this version number.
 
 @returns Stimulus version.
 */
-(NSUInteger)version;

@end

/**
AUIContinuousStimulus defines an added method for stimuli that allow
 continuous stimulus generation (i.e. not in 'blocks').
 */
@protocol AUIContinuousStimulus<AUIStimulus>
/*
 If implemented, stimulus provides support for continuous stimulus generation.
 @returns Next stimulus sample
 */
-(double)nextSample;
@end