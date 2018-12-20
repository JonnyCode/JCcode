function [ReturnConditionA, ReturnConditionB] = SplitCondition(CellInfo, Condition, Fraction)

NumEpochs = length(Condition.ExcludeEpochs);	

EndEpoch = floor(NumEpochs * Fraction);
ReturnConditionA.EpochNumbers = Condition.EpochNumbers(1:EndEpoch);
ReturnConditionB.EpochNumbers = Condition.EpochNumbers(EndEpoch + 1:NumEpochs);
ReturnConditionA.ExcludeEpochs = Condition.ExcludeEpochs(1:EndEpoch);
ReturnConditionB.ExcludeEpochs = Condition.ExcludeEpochs(EndEpoch + 1:NumEpochs);
ReturnConditionA.UserInfo.StimulusStartTime = Condition.UserInfo.StimulusStartTime(1:EndEpoch);
ReturnConditionB.UserInfo.StimulusStartTime = Condition.UserInfo.StimulusStartTime(EndEpoch + 1:NumEpochs);
ReturnConditionA.EpochTime = Condition.EpochTime(1:EndEpoch);
ReturnConditionB.EpochTime = Condition.EpochTime(EndEpoch + 1:NumEpochs);
ReturnConditionA.EpochData.Offset = Condition.EpochData.Offset(1:EndEpoch);
ReturnConditionB.EpochData.Offset = Condition.EpochData.Offset(EndEpoch + 1:NumEpochs);
ReturnConditionA.SearchCrit = Condition.SearchCrit;
ReturnConditionA.SearchPara = Condition.SearchPara;
ReturnConditionB.SearchCrit = Condition.SearchCrit;
ReturnConditionB.SearchPara = Condition.SearchPara;

for epoch = 1:EndEpoch
	ReturnConditionA.EpochData.Data(epoch, :) = Condition.EpochData.Data(epoch, :);
end

for epoch = EndEpoch+1:NumEpochs
	ReturnConditionB.EpochData.Data(epoch - EndEpoch, :) = Condition.EpochData.Data(epoch, :);
end

ReturnConditionA = AverageAndVariance(CellInfo, ReturnConditionA);
ReturnConditionB = AverageAndVariance(CellInfo, ReturnConditionB);
