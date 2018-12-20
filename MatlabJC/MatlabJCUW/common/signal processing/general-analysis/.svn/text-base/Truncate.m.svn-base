function ReturnCondition = Truncate(Condition, RespLength)
% Truncate
%
% Truncate data in Condition at length specified. 

NumEpochs = length(Condition.ExcludeEpochs);	
ReturnCondition = Condition;
ReturnCondition = rmfield(ReturnCondition, 'EpochData');
ReturnCondition.EpochPts = RespLength;
EpochPts = Condition.EpochPts;
for epoch = 1:NumEpochs
	temp = Condition.EpochData(epoch, 1:RespLength);
	ReturnCondition.EpochData(epoch, :) = temp;
end

ReturnCondition = AverageAndVariance(ReturnCondition);
