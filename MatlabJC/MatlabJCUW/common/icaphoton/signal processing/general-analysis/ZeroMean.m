function [ReturnedCellInfo, ReturnedCondition] = ZeroMean(CellInfo, Condition, StartPt, Points)
%
% [ReturnedCellInfo, ReturnedCondition] = ZeroMean(CellInfo, Condition)
%
%  Takes data and subtracts any offset.
%  Created:  FMR  ?
%	Revised:   GDF  9/14/01
%				-made compatible with new data format
%			FMR 9/8/03
%				- added option to use only subset of responses
				
ReturnedCellInfo = CellInfo;
NumEpochs = length(Condition.ExcludeEpochs);	
ReturnedCondition = Condition;
[EpochData, Offset] = GetEpochData(CellInfo, Condition);
EpochPts = size(EpochData, 2);
FieldFlag = isfield(Condition, 'EpochData');
TempData = EpochData(:, StartPt:StartPt+Points-1);
AveRow = mean(TempData, 2);
for epoch = 1:NumEpochs
	EpochData(epoch,:) = EpochData(epoch, :) - AveRow(epoch);
end	

clear ReturnedCondition.EpochData.Data

if (FieldFlag == 1)
	ReturnedCondition.EpochData.Data = EpochData;
end

if (FieldFlag == 0)
	for cnt = 1:NumEpochs
		CurrentEpoch = Condition.EpochNumbers(cnt) + 1;
		CellInfo.EpochData.Data{CurrentEpoch} = EpochData(cnt);
	end
end

ReturnedCondition = AverageAndVariance(CellInfo, ReturnedCondition);
