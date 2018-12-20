function ReturnedCondition = AddNoiseToCondition(CellInfo, Condition, NoiseAmp);
% AddNoiseToCondition
%
%	Adds zero-mean gaussian white noise to each epoch in condition.  The standard
% 	deviation of the noise is specified by NoiseAmp.


EpochPts = length(Condition(1).AverageResponse);
NumConds = length(Condition);

for cond = 1:NumConds
	Condition(cond).EpochData.Data = Condition(cond).EpochData.Data + random('norm', 0, NoiseAmp, size(Condition(cond).EpochData.Data, 1), size(Condition(cond).EpochData.Data, 2));
end

ReturnedCondition = Condition;
