function Projections = Calculate2AFCResponseProjection3(CellInfo, ConditionA, ConditionB, NoiseCondition, cond, FreqCutoff, RespLength, Shuffles)
%
%	 Projections = Calculate2AFCResponseProjection2(CellInfo, ConditionA, ConditionB, NoiseCondition, cond, FreqCutoff, RespLength, Shuffles)
%
% CalculateResponseProjection3
%
%	Compute  response projections along discriminant for two 
% 	sets of data in 2IFC experiment. 
%
%	11/03 FMR Created


NumIterations = floor(length(ConditionA(cond).ExcludeEpochs)/2);	% number of data epochs in train and test

PCorrect = 0;
EpochList = ConditionA(cond).EpochNumbers;

GoodValues = 0;
NoiseValues = 0;
	
% find discriminant from training data
Discriminant = FindDiscriminant2(CellInfo, ConditionA, ConditionB, cond, FreqCutoff, RespLength, EpochList);
NormFact = sum(Discriminant .* Discriminant);
		
[respA, Offset] = GetGoodEpochData(CellInfo, ConditionA(cond));
[respB, Offset] = GetGoodEpochData(CellInfo, ConditionB(cond));
[noise, Offset] = GetGoodEpochData(CellInfo, NoiseCondition);

% apply to each of test data epochs
clear Projections;
Projections.RespAProjection = respA * Discriminant' / NormFact;	% response projections
Projections.RespBProjection = respB * Discriminant' / NormFact; 	% response projections for second data set
Projections.NoiseProjection = noise * Discriminant' / NormFact; 	% response projections for second data set

