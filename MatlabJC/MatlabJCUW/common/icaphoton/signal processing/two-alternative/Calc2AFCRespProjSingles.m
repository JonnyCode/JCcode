function Projections = Calc2AFCRespProjSingles(CellInfo, ConditionA, NoiseCondition, cond, FreqCutoff, RespLength, Shuffles)
%
%	 Projections = Cal2AFCRespProjSingles(CellInfo, ConditionA, NoiseCondition, cond, FreqCutoff, RespLength, Shuffles)
%
% CalculateResponseProjection2
%
%	Compute  response projections along discriminant for two 
% 	sets of data in 2IFC experiment.  First split data into training and testing sets.
%	Then calculate discriminant for training data.  Then apply discriminant
%	to test data and compute projections.
%
%	4/30/00 FMR Created
%	10/9/00 FMR Revised to handle Condition structures
%	10/8/01 GDF Revised to handle new data structures
%	10/01	GDF Revised to extract single photon responses

NumIterations = floor(length(ConditionA(cond).ExcludeEpochs)/2);	% number of data epochs in train and test
RespAProjection(1:Shuffles * NumIterations) = 0;	% response projections
RespBProjection(1:Shuffles * NumIterations) = 0; 	% response projections for second data set
NoiseProjection(1:Shuffles * NumIterations) = 0; 	% response projections for second data set

PCorrect = 0;
EpochList = 1:NumIterations*2;

GoodValues = 0;
NoiseValues = 0;
for Shuffle = 0:Shuffles-1

	% generate train and test data sets	
	ShuffledList = ShuffleList(EpochList);
	Training = randperm(NumIterations);
	Test = ShuffleList(NumIterations+1:2*NumIterations);
	
	% find discriminant from training data
	Discriminant = FindDiscrimForSingles(CellInfo, ConditionA, NoiseCondition, cond, FreqCutoff, RespLength, Training);
	NormFact = sum(Discriminant .* Discriminant);
		
	% apply to each of test data epochs
	for cnt = 1:NumIterations
		if (ConditionA(cond).ExcludeEpochs(Test(cnt)) == 0)
			if cnt > length(NoiseCondition.ExcludeEpochs)
				noisecnt = cnt - length(NoiseCondition.ExcludeEpochs);
			else
				noisecnt = cnt;
			end
			if Test(noisecnt) > length(NoiseCondition.ExcludeEpochs)
				noiseTest = Test(noisecnt) - length(NoiseCondition.ExcludeEpochs);
			else
				noiseTest = Test(noisecnt);
			end			
			if (NoiseCondition.ExcludeEpochs(noiseTest) == 0)
				GoodValues = GoodValues + 1;
				respA = ConditionA(cond).EpochData.Data(Test(cnt),:);
				respB = NoiseCondition.EpochData.Data(noiseTest,:);
				% compute projections of responses onto discriminant and update probability correct
				RespAVal = sum(respA .* Discriminant) / NormFact;
				RespBVal = sum(respB .* Discriminant) / NormFact;
				RespAProjection(GoodValues) = RespAVal;
				RespBProjection(GoodValues) = RespBVal;
			end
		end
	end

	% calculation noise projections
	for cnt = 1:length(NoiseCondition.ExcludeEpochs)
		if cnt > length(NoiseCondition.ExcludeEpochs)
			noisecnt = cnt - length(NoiseCondition.ExcludeEpochs);
		else
			noisecnt = cnt;
		end
		if (NoiseCondition.ExcludeEpochs(noisecnt) == 0)
			NoiseValues = NoiseValues + 1;
			noiseresp = NoiseCondition.EpochData.Data(noisecnt,:);		
			NoiseVal = sum(noiseresp .* Discriminant) / NormFact;
			NoiseProjection(NoiseValues) = NoiseVal;
		end
	end
end

	
clear Projections;
Projections.RespAProjection(1:GoodValues) = RespAProjection(1:GoodValues);	% response projections
Projections.RespBProjection(1:GoodValues) = RespBProjection(1:GoodValues); 	% response projections for second data set
Projections.NoiseProjection(1:NoiseValues) = NoiseProjection(1:NoiseValues); 	% response projections for second data set

