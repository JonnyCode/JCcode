function ReturnCondition = AverageAndVariance(ConditionA, ConditionB)

NumEpochs = length(ConditionA.ExcludeEpochs);	
ReturnCondition = ConditionA;

for epoch = 1:NumEpochs
	if (ConditionB.ExcludeEpochs(epoch) == 1)
		ReturnCondition.ExcludeEpochs(epoch) == 1;
	end
	ReturnCondition.EpochData.Data(epoch, :) = ConditionA.EpochData.Data(epoch, :) - ConditionB.EpochData.Data(epoch, :);
end

