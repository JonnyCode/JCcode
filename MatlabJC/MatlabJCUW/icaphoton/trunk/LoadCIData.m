function ReturnedCellInfo = LoadCIData(CellInfo,varargin)
% ReturnedCellInfo = LoadCIData(CellInfo,condition)
% version 2.1
% condition is optional - 1 is epochcondition, 2 is familycondition, 3 for both
% - does so without any interaction, if conition is empty, default is
% interactive.
%
% version 2 combines old loadcidata and decimate
% Function to load epoch data into the structure CellInfo.  Structure must 
% be loaded in the workspace.  ReturnedCellInfo can also be called CellInfo. Works 
% for both Epoch and Family Conditions, and asks the user about loading the data 
% before going through either one. The 
% program gets the epoch numbers for each condition, zeros the data, loads the 
% offset into CellInfo.EpochData.Offset, and decimates the data by the number
% in the field CellInfo.(x)Condition.Decimate for each Condition.  Decimate 
% reshapes the data into columns and takes the mean, so that no data is actually
% lost.  In other words if you decimate by two, then we take the average of every
% two points for the new data structure. The data is then loaded into 
% ReturnedCellInfo.EpochData.Data.  To 
% remove the field use the script RmCIData.  To load data in files that have more 
% than one segment, use LoadSCIData.

% Created as loadSCData February, 2001 MKMK
% new version, LoadCIData released August, 2001 MKMK
% version 2, LoadCIData released February, 2002 MKMK
% Changed to a function, uses CellInfo from the workspace instead of a file
% version 2.1, added input condition, so function can be run w/o any interaction
% January, 2003 MKMK
% version 2.2, added a check for the path.
% April, 2005 MKMK removed insult, only give warning about data already
% being there once, and don't bother zeroing, decimating, and putting it there again.

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
[fileptr, ecode] = ITCInitializeAnalysis(2e6, fp);
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
dataCheck = 0;
if Familyes == 1
	Numconditions = length(CellInfo.FamilyCondition);
	for i = 1:Numconditions
		Epochnums = CellInfo.FamilyCondition(i).EpochNumbers;
		Numepochs = length(Epochnums);
		decimatefactor = CellInfo.FamilyCondition(i).DecimatePts;
        for j = 1:Numepochs
            Epochnum = Epochnums(j);
            % check to see if there is already data, if so, skip this epoch
            if isempty(CellInfo.EpochData.Data{Epochnum + 1})
                [epochdata,ecode] = ITCReadEpoch(Epochnum, segnum, fileptr);
                if ecode
                    fprintf(1, 'error1 reading from file \n');
                    error
                end
                [prepoints,ErrorFlag] = ITCGetStmPrePts(Epochnum, segnum, 0, fileptr);
                if ErrorFlag
                    fprintf(1, 'error1 reading from file \n');
                    error
                end
                % zero the data
                trial = epochdata(1:prepoints + 1);
                baselineaverage = mean(trial);
                epochdata = epochdata - baselineaverage;
                epochlength = length(epochdata);
                % decimate the data
                otherfactor = epochlength/decimatefactor;
                newshape = reshape(epochdata,decimatefactor,otherfactor);
                newdata = mean(newshape,1);
                % read the new data in, and record the offset
                CellInfo.EpochData.Data{Epochnum + 1} = newdata;
                CellInfo.EpochData.Offset(Epochnum + 1) = baselineaverage;
            elseif dataCheck == 0
                fprintf(1, 'There is already data here, skipping epoch \n')
                dataCheck = 1;
            end
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
	for i = 1:Numconditions
		Epochnums = CellInfo.EpochCondition(i).EpochNumbers;
		Numepochs = length(Epochnums);
		decimatefactor = CellInfo.EpochCondition(i).DecimatePts;
        for j = 1:Numepochs
            Epochnum = Epochnums(j);
            % check to see if there is already data, if so, skip this epoch
            if isempty(CellInfo.EpochData.Data{Epochnum + 1})
                [epochdata,ecode] = ITCReadEpoch(Epochnum, segnum, fileptr);
                if ecode
                    fprintf(1, 'error1 reading from file \n');
                    error
                end
                [prepoints,ErrorFlag] = ITCGetStmPrePts(Epochnum, segnum, 0, fileptr);
                if ErrorFlag
                    fprintf(1, 'error1 reading from file \n');
                    error
                end
                % zero the data 
                trial = epochdata(1:prepoints + 1);
                baselineaverage = mean(trial);
                epochdata = epochdata - baselineaverage;
                epochlength = length(epochdata);
                % decimate the data
                otherfactor = epochlength/decimatefactor;
                newshape = reshape(epochdata,decimatefactor,otherfactor);
                newdata = mean(newshape,1);
                % read the new data in, and record the offset
                CellInfo.EpochData.Data{Epochnum + 1} = newdata;
                CellInfo.EpochData.Offset(Epochnum + 1) = baselineaverage;
            elseif dataCheck == 0
                fprintf(1, 'There is already data here, skipping epoch \n')
                dataCheck = 1;
            end
        end
	end
end
ReturnedCellInfo = CellInfo;
