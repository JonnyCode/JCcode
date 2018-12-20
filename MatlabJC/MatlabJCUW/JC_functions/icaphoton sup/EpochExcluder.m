function CellInfo = EpochExcluder(CellInfo,epochCond_num,epochsUwant) ; 

% This function will make sure the epochs you want are in CellInfo and
% exclude the epochs you don't want
% Input = CellInfo file, epochcondtion number, vector of epochs you
% want

% JC 9/4/07

% TO BE USED AS A SCRIPT
% epochCond_num = ; enter the epochcondition number
% epochsUwant = [] ; % enter the epochs you want to analyze

[a,b,c] = intersect(CellInfo.EpochCondition(epochCond_num).EpochNumbers,epochsUwant) ; % a= the intersection, b=index of CellInfo... that are in the intesrsection, c=index of epochs that are in the intersection

if length(epochsUwant)~=length(a) ; % if any one of the epochs you wnat is not in the CellInfo file.... 
    disp('CellInfo is missing one or more of the epochs you wanted, apologies all around') ; % display this message
else 
    CellInfo.EpochCondition(epochCond_num).ExcludeEpochs = ones(1,length(CellInfo.EpochCondition(epochCond_num).ExcludeEpochs)) ; % exculde all epochs
    CellInfo.EpochCondition(epochCond_num).ExcludeEpochs(b) = 0 ;  % don't exclude the epochs you wanted
end
