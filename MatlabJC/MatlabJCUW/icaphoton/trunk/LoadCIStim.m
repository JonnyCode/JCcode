function ReturnedCellInfo = LoadCIStim(CellInfo,varargin)
% ReturnedCellInfo = LoadCIStim(CellInfo,condition)
% version 1
% condition is optional - 1 is epochcondition, 2 is familycondition, 3 for
% both - does so without any interaction
%
% Function to load epoch stimulus into the structure CellInfo.  Structure
% must be loaded in the workspace.  ReturnedCellInfo can also be called
% CellInfo. Works for both Epoch and Family Conditions, and asks the user
% about loading the data before going through either one. The program gets
% the epoch numbers for each condition, and decimates the data by the number in
% the field CellInfo.(x)Condition.Decimate for each Condition.  Decimate
% reshapes the data into columns and takes the mean, so that no data is
% actually lost.  In other words if you decimate by two, then we take the
% average of every two points for the new data structure. The data is then
% loaded into ReturnedCellInfo.EpochData.Data.
%
%
% December, 2004 MKMK

varargin{:};
if nargin<2
	condition = [];
else
	condition = varargin{1};
end
fp = CellInfo.CellFile;
if ~isempty(findstr(':',fp))
    if ~strcmp(version,'5.2.1.1421')
        disp(fp)
        fp = input('What is the location plus name of the CONVERTED (for OS X) data file? ','s');
    end
end
[fileptr, ecode] = ITCInitializeAnalysis(500000, fp);
[NumEpochs, ecode5] = ITCGetNumberEpochs(fileptr);
if ecode5~=0
	error
end
CellInfo.EpochData.Data = cell(NumEpochs,1);
CellInfo.EpochData.Offset = zeros(NumEpochs,1);
segnum = 0;
if isempty(condition)
	fprintf(1, 'Default is yes \n')
	Familycheck = input('Go through FamilyConditions? y or n ','s');
	if isempty(Familycheck)
		Familyes = 1;
	else
		Familyes = strcmp(Familycheck,'y');
	end
else
	if condition == 1
		Epochyes = 1;
		Familyes = 0;
	elseif condition == 2
		Familyes = 1;
		Epochyes = 0;
	elseif condition == 3
		Familyes = 1;
		Epochyes = 1;
	else
		error('Please use 1,2 or 3 for the condition input, see help for further info')
	end
end	
if Familyes == 1
	Numconditions = length(CellInfo.FamilyCondition);
	for i = 1:Numconditions
		Epochnums = CellInfo.FamilyCondition(i).EpochNumbers;
		Numepochs = length(Epochnums);
		decimatefactor = CellInfo.FamilyCondition(i).DecimatePts;
		for j = 1:Numepochs
			Epochnum = Epochnums(j);
            [StimData, ecode] = ITCReadEpochStm(Epochnum, segnum, fileptr);
			epochlength = length(StimData);
			otherfactor = epochlength/decimatefactor;
			newshape = reshape(StimData,decimatefactor,otherfactor);	
			newdata = mean(newshape,1);
			CellInfo.EpochData.StimData{Epochnum + 1} = newdata;
		end
	end
end
if isempty(condition)
	fprintf(1, 'Default is yes \n')
	Epochcheck = input('Go through EpochConditions? y or n ','s');
	if isempty(Epochcheck)
		Epochyes = 1;
	else
		Epochyes = strcmp(Epochcheck,'y');
	end
end
if Epochyes == 1
	Numconditions = length(CellInfo.EpochCondition);
	maxEpoch = 0;
	for i = 1:Numconditions
		maxEpoch = max([maxEpoch max(CellInfo.EpochCondition(i).EpochNumbers)]);
	end
	
	for i = 1:Numconditions
		Epochnums = CellInfo.EpochCondition(i).EpochNumbers;
		Numepochs = length(Epochnums);
		decimatefactor = CellInfo.EpochCondition(i).DecimatePts;
		for j = 1:Numepochs
			Epochnum = Epochnums(j);
			[StimData, ecode] = ITCReadEpochStm(Epochnum, segnum, fileptr);
			epochlength = length(StimData);
			otherfactor = epochlength/decimatefactor;
			newshape = reshape(StimData,decimatefactor,otherfactor);	
			newdata = mean(newshape,1);
			if(~isfield(CellInfo.EpochData, 'StimData'))
				CellInfo.EpochData.StimData=cell(1,maxEpoch+1);
			end
			if isempty(CellInfo.EpochData.StimData{Epochnum + 1})
				CellInfo.EpochData.StimData{Epochnum + 1} = newdata;
			else
				fprintf(1,'data there already! ... overwriting');
				CellInfo.EpochData.StimData{Epochnum + 1} = newdata;
			end
		end
	end
end
ReturnedCellInfo = CellInfo;
