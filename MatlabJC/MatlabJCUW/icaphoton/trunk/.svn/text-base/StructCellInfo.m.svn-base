%	CellInfo structure helpfile 
%
%	CellInfo.CellFile 				File Name of original data file, incl. path
%	CellInfo.Comment 				Filled in by user
%	CellInfo.Label 					Filled in by user
%	CellInfo.CellType 				Filled in by user
%	CellInfo.Rig 					Filled in by user
%	CellInfo.OutputScaleFactor		How to scale the output for this setting
%	CellInfo.OutputConfiguration    Which output corresponds to which output setting
%	CellInfo.NDFConfiguration		What NDF configuration was used for this setting
%	CellInfo.UserInfo				Can be structure, user defined
%	CellInfo.EpochCondition			Structure, defined by the epochs all
%										exactly the same as defined by user
%	CellInfo.FamilyCondition		Structure, defined by groups of epochs
%										that increase in value of one variable
%	CellInfo.EpochData				Cell Structure of all the data in the file
%										Use LoadCIData and RmCIData to load/remove 
%
%  		These I am calling Generic Condition because they can be found in both the 
% 		EpochCondition and the FamilyCondition structures within CellInfo:
%
%	GenericCondition.Segnum			Segment number analyzed	    
%	GenericCondition.Epochnumbers  	List of epoch numbers as selected in Epochs
%  	GenericCondition.ExcludeEpochs	List of 0's and 1's, length of Epochnumbers, 
%										0's deselected during analysis
%	GenericCondition.Comment		Filled in by user
%	GenericCondition.Plotpref		plotpref is defined in icaphoton, user 	
%										sets preference for plotting 
%	GenericCondition.Searchcrit		Criteria used to choose the epochs (ie. same values 
%										for all epochs)								
%	GenericCondition.Searchpara		parameters (actual values of criteria) used 
%										to select epochs
%   GenericCondition.UserInfo 		Can be structure, user defined.  If use LEDCalibration
%										routine, one field will be StimulusAmp, the amp 
%										converted for that particular output channel.  For 
%										FamilyCondition this will also have a new version of 
%										familystep
%   GenericCondition.DecimatePts    Number that represents how much the data is decimated by.
%										Data is decimated automatically when data is loaded 
%										into the structure using LoadCIData or LoadSCIData.
%										See help LoadCIData or LoadSCIData for more info.
%
%		These are unique to FamilyCondition
%
%	FamilyCondition.FamilyFlag		size of EpochNumbers, gives each family a sequential number									
%	FamilyCondition.FamilyStep		values of step parameter for epochs in EpochNumbers
%	FamilyCondition.FamilyStepGuide	name of parameter that steps (ie. amp)
%
%		These can be empty if all values of the epochs besides the step are the same:
%	FamilyCondition.FamilyCues		values of cue parameter for epochs in EpochNumbers
% 	FamilyCondition.FamilyCueGuide	name of parameter(s) that change from family to family
%
%
