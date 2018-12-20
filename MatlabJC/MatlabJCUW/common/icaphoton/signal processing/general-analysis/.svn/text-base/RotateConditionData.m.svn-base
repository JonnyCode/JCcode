function ReturnCondition = RotateConditionData(Condition, RotatePoints)
%
%	ReturnCondition = RotateResponseVectors(Condition, RotatePoints)
%
%	Takes a condition and rotates the data by RotatePoints

NumEpochs = length(Condition.EpochNumbers);		% Get the number of epochs
ReturnCondition = Condition;

% go through each epoch and extract appropriate A and B vectors
for cnt = 1:NumEpochs
	TempResp = Condition.EpochData.Data(cnt, :);
	TempResp = RotateVector(RotatePoints, TempResp);
	ReturnCondition.EpochData.Data(cnt,:) = TempResp;
end
