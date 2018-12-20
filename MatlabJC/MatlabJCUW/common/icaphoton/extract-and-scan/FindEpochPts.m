function EpochPts = FindEpochPts(CellInfo, Condition)
%  EpochPts = FindEpochPts(CellInfo, Condition)
%
%  	This function returns the number of data points
%	in a given condition (either EpochCondition or FamilyCondition).
%	Inputs are a CellInfo and a Condition.
%  	It will get the number of epoch points from a CellInfo
%  	if there is no data in the Condition.
%	
%	Created ? GDF

FieldFlag = isfield(Condition, 'EpochData');
if FieldFlag == 0
	% Get the number of EpochPnts from the EpochData in CellInfo
	CurrentEpoch = Condition.EpochNumbers(1) + 1;
	EpochPts = length(CellInfo.EpochData.Data{CurrentEpoch});
end

if FieldFlag == 1
	% Get the number of EpochPnts from the EpochData in the EpochCondition
	EpochPts = length(Condition.EpochData.Data(1,:));
end	
