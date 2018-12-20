% Do2AFCDiscrimination
%
%	 PCorrect = Do2AFCDiscrimination(CellInfo, RespConditionA, RespConditionB, cond, shuffles, FreqCutoff, RespLength, Verbose)
%
%	Do two alternative forced choice discrimination between A and B vectors.
%	Repeat for several random shufflings (pairings) of A and B responses.
%
% Created 10/00 FMR
%
% Revisions
%	10/22/00 FMR removed mean correction
%	7/15/01  FMR created from two-interval discrimination
%	9/14/01  GDF modified to handle new data format

function PCorrect = Do2AFCDiscrimination(CellInfo, RespConditionA, RespConditionB, cond, shuffles, FreqCutoff, RespLength, Verbose)

PCorrect = 0;						% probability correct
NumIterations = floor(length(RespConditionA(cond).ExcludeEpochs)/2);
list = 1:NumIterations*2;					% list for resampling
NumComparisons = 0;
EpochPts = length(RespConditionA(cond).AverageResponse);
FlashPtA = round(RespConditionA(cond).UserInfo.StimulusStartTime(1) / (FindSearchPara(RespConditionA(cond), 'SampInterv') * 1e-6 * RespConditionA(cond).DecimatePts));
FlashPtB = round(RespConditionB(cond).UserInfo.StimulusStartTime(1) / (FindSearchPara(RespConditionB(cond), 'SampInterv') * 1e-6 * RespConditionB(cond).DecimatePts));

% resample data shuffles times, using half as learn and half as test.
for shuffle = 1:shuffles
	
	if (Verbose)
		fprintf(1, 'Shuffle %d\n', shuffle);
	end
	
	% shuffle list and make test and learn data sets
	ShuffledList = randperm(NumIterations * 2);
	LearnList = randperm(NumIterations);
	TestList = ShuffleList(NumIterations+1:2*NumIterations);
		
	% compute discriminant from first half of list
	Discriminant = FindDiscriminant(CellInfo, RespConditionA, RespConditionB, cond, FreqCutoff, RespLength, LearnList);
	plot(Discriminant);
	
	% apply discriminant to test list -- first do A vectors in first half list
	for cnt = 1:floor(NumIterations/2)
		if (RespConditionA(cond).ExcludeEpochs(TestList(cnt)) == 0)
			% compute projections of responses onto discriminant and update probability correct
			resp = RespConditionA(cond).EpochData.Data(TestList(cnt), :);
			if (RespLength < EpochPts)
				resp(1:FlashPtA) = 0;
				resp(FlashPtA + RespLength:EpochPts) = 0;
			end
			RespVal = sum(resp .* Discriminant);
			if (RespVal > 0)
				PCorrect = PCorrect + 1;
			end
			if (RespVal == 0)
				PCorrect = PCorrect + 0.5;
			end
			NumComparisons = NumComparisons + 1;
		end
		if (Verbose)
			if (rem(cnt, 10) == 0)
				fprintf(1, '\t%d: Probability Correct = %d\n', cnt, PCorrect/NumComparisons);
			end
		end	
	end
	% now do B vectors in second half of test list
	for cnt = floor(NumIterations/2)+1:NumIterations
		if (RespConditionB(cond).ExcludeEpochs(TestList(cnt)) == 0)
			% compute projections of responses onto discriminant and update probability correct
			resp = RespConditionB(cond).EpochData.Data(TestList(cnt), :);
			if (RespLength < EpochPts)
				resp(1:FlashPtB) = 0;
				resp(FlashPtB + RespLength:EpochPts) = 0;
			end
			RespVal = sum(resp .* Discriminant);
			if (RespVal < 0)
				PCorrect = PCorrect + 1;
			end
			if (RespVal == 0)
				PCorrect = PCorrect + 0.5;
			end
			NumComparisons = NumComparisons + 1;
		end
		if (Verbose)
			if (rem(cnt, 10) == 0)
				fprintf(1, '\t%d: Probability Correct = %d\n', cnt, PCorrect/NumComparisons);
			end
		end	
	end
end

PCorrect = PCorrect / NumComparisons;
fprintf(1, 'Probability Correct = %d\n', PCorrect);
