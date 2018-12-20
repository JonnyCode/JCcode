function ReturnCondition = InterpolateCondition(Condition, InterpFactor);
%
%	function ReturnCondition = InterpolateCondition(Condition, InterpFactor);
%
%	Interpolate condition expands the data in a condition.  It is the opposite of
%	decimate.  InterpFactor is the increase in sampling rate - i.e. a factor
% 	of 2 will result in twice as many data points.

ReturnCondition = Condition;
EpochData = Condition.EpochData.Data;
[NumEpochs, EpochPts] = size(EpochData);
for epoch = 1:NumEpochs
	clear temp newtemp
	temp = EpochData(epoch,:);
	newtemp = interp(temp, InterpFactor);
	NewEpochData(epoch,:) = newtemp;
end
ReturnCondition.EpochData = rmfield(ReturnCondition.EpochData,'Data');
ReturnCondition.EpochData.Data = NewEpochData;
ReturnCondition = rmfield(ReturnCondition,'AverageResponse');
ReturnCondition = rmfield(ReturnCondition,'VarianceResponse');
ReturnCondition.AverageResponse = mean(NewEpochData);
ReturnCondition.VarianceResponse = var(NewEpochData);
ReturnCondition.DecimatePts = ReturnCondition.DecimatePts / InterpFactor;
