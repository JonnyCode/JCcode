function [ReturnCellInfo, ReturnCondition] = ZeroBaseLine(CellInfo, Condition, BaseLinePts)
%
%	function ReturnCondition = ZeroBaseLine(CellInfo, Condition, BaseLinePts)
%
%	This function averages over the baseline values and then subtracts this
%	average from the rest of the data
%	Created GDF 12/01

ReturnCondition = Condition;
ReturnCellInfo = CellInfo;
if isfield(Condition, 'EpochData')
	[NumEpochs EpochPts] = size(Condition.EpochData.Data);
	TempData(:,:) = Condition.EpochData.Data(:,1:BaseLinePts);
	AveBaseLines = mean(TempData, 2);
	for epoch = 1:NumEpochs;
		NewData(epoch,:) = ReturnCondition.EpochData.Data(epoch,:)-AveBaseLines(epoch);
	end
	ReturnCondition.EpochData = rmfield(ReturnCondition.EpochData, 'Data');
	ReturnCondition.EpochData.Data = NewData;
else
	fprintf(1,'You are screwed because this function has not been completed')
end
