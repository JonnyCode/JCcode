function [ReturnConditionA, ReturnConditionB] = ShuffleEpochs(ConditionA, ConditionB)
% ShuffleEpochs
%
% Randomly shuffle epochs in two conditions.

NumEpochs = length(ConditionA.ExcludeEpochs);	
list = 1:NumEpochs;
ShuffledList = ShuffleList(list);

ReturnConditionA = ConditionA;
ReturnConditionB = ConditionB;

for epoch = 1:NumEpochs
	temp = ConditionA.EpochData.Data(ShuffledList(epoch), :);
	ReturnConditionA.EpochData.Data(epoch, :) = temp;
	temp = ConditionB.EpochData.Data(ShuffledList(epoch), :);
	ReturnConditionB.EpochData.Data(epoch, :) = temp;
	ReturnConditionA.ExcludeEpochs(epoch) = ConditionA.ExcludeEpochs(ShuffledList(epoch));
	ReturnConditionB.ExcludeEpochs(epoch) = ConditionB.ExcludeEpochs(ShuffledList(epoch));
	ReturnConditionA.EpochNumbers(epoch) = ConditionA.EpochNumbers(ShuffledList(epoch));
	ReturnConditionB.EpochNumbers(epoch) = ConditionB.EpochNumbers(ShuffledList(epoch));
end
