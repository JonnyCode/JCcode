function ReturnCondition = ApplyNonlinearity(CellInfo, Condition, Mean, SD)
% ApplyNonlinearity
%
% 	Apply multiplicative nonlinearity to all data on condition.  Nonlinearity
%	is cumulative gaussian, and each point in each epoch is multiplied by 
%	corresponding weight of the cumulative gaussian.
%

NumEpochs = length(Condition.EpochNumbers);	
ReturnCondition = Condition;
EpochPts = length(Condition.AverageResponse);

for epoch = 1:NumEpochs
	temp(1:EpochPts) = Condition.EpochData.Data(epoch, 1:EpochPts);
	ReturnCondition.EpochData.Data(epoch, 1:EpochPts) = normcdf(temp, Mean, SD) .* temp;
end

ReturnCondition = AverageAndVariance(CellInfo, ReturnCondition);
