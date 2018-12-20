function ReturnCondition = DecimateCondition(CellInfo, Condition, DecimatePoints)
% DecimateCondition
%
%	Decimate data in Condition.

ReturnCondition = Condition;
ReturnCondition.EpochData = rmfield(Condition.EpochData, 'Data');
EpochData = Condition.EpochData.Data;
for epoch = 1:size(EpochData, 1)
	 NewEpochData(epoch, :) = decimate(EpochData(epoch, :), DecimatePoints);
end
ReturnCondition.EpochData.Data = NewEpochData;
ReturnCondition = AverageAndVariance(CellInfo, ReturnCondition);
ReturnCondition.DecimatePts = DecimatePoints;
