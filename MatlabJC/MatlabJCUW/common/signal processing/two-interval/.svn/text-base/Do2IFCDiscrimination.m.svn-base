% Do2IFCDiscrimination
%	Do two interval forced choice discrimination between A and B vectors.
%	Repeat for several random shufflings (pairings) of A and B responses.
%
% Created 10/00 FMR
%
% Revisions
%	10/22/00 FMR removed mean correction
%	10/09/02 FMR updated to new data format

function PCorrect = DoDiscrimination(CellInfo, RespConditionA, RespConditionB, cond, shuffles, FreqCutoff, RespLength, Verbose)

PCorrect = 0;						% probability correct
NumIterations = floor(length(RespConditionA(cond).ExcludeEpochs)/2);
list = 1:NumIterations*2;					% list for resampling

NumComparisons = 0;

% resample data shuffles times, using half as learn and half as test.
for shuffle = 1:shuffles
	
	if (Verbose)
		fprintf(1, 'Shuffle %d\n', shuffle);
	end
	
	% shuffle list and make test and learn data sets
	ShuffledList = ShuffleList(list);
	LearnList = ShuffledList(1:NumIterations);
	TestList = ShuffledList(NumIterations+1:2*NumIterations);
	RotatedTestList = RotateVector(1, TestList);
	
	% compute discriminant from first half of list
	Discriminant = FindDiscriminant(CellInfo, RespConditionA, RespConditionB, cond, FreqCutoff, RespLength, LearnList);
	plot(Discriminant);
	
	% apply discriminant to second half	
	for cnt = 1:NumIterations
		if (RespConditionA(cond).ExcludeEpochs(TestList(cnt)) == 0)
			if (RespConditionB(cond).ExcludeEpochs(TestList(cnt)) == 0)
				% compute projections of responses onto discriminant and update probability correct
				respA = RespConditionA(cond).EpochData.Data(TestList(cnt), :);
				respB = RespConditionB(cond).EpochData.Data(RotatedTestList(cnt), :);
				RespAVal = sum(respA .* Discriminant);
				RespBVal = sum(respB .* Discriminant);
				if (RespAVal > RespBVal)
					PCorrect = PCorrect + 1;
				end
				if (RespAVal == RespBVal)
					PCorrect = PCorrect + 0.5;
				end
				NumComparisons = NumComparisons + 1;
			end
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
