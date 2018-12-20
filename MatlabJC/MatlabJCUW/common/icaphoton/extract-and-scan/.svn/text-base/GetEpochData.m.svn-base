function [EpochData, Offset] = GetEpochData(CellInfo, Condition)
%  GetEpochData.m
%
%  [EpochData, Offset] = GetEpochData(CellInfo, Condition)
%
%  Gets the EpochData from either a CellInfo or from a Condition
%  GetEpochData outputs all the data at once in a matrix.  
%  All analysis on the data should be geared to analyzing the data in this
%  manner, rather than one vector at a time.

FieldFlag = isfield(Condition, 'EpochData');
NumEpochs = length(Condition.ExcludeEpochs);
if FieldFlag == 0
	EpochNumbers = Condition.EpochNumbers + 1;
	for epoch = 1:NumEpochs
		EpochData(epoch,:) = CellInfo.EpochData.Data{EpochNumbers(epoch), Condition.SegNum+1};
		Offset(epoch) = CellInfo.EpochData.Offset(EpochNumbers(epoch),Condition.SegNum+1);
	end
end
if FieldFlag == 1
	EpochData = Condition.EpochData.Data;
	Offset = Condition.EpochData.Offset;
end
