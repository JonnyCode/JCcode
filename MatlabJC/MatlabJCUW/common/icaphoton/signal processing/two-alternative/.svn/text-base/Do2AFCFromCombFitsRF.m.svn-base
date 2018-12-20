function DiscrimParams = Do2AFCFromCombFitsRF(CombParameters, RFWeight, NumTrials, MeanPhotoisoms)

% create two sets of samples from comb distribution --- for flash conditions A and B --- 
% and do discrimination to compare to true measured responses
NumTestConds = length(MeanPhotoisoms);
clear LinPCorrect;
clear NonLinPCorrect;
LinPCorrect(1:NumTestConds) = 0;
NonLinPCorrect(1:NumTestConds) = 0;
clear respA;
clear respB;

for cond=1:NumTestConds 
	CombParameters.RateA = MeanPhotoisoms(cond);
	CombParameters.RateB = MeanPhotoisoms(cond);
	for resp=1:NumTrials;
		[respA, respB] = SampleFromCombDistribution(CombParameters, length(RFWeight)^2, 1);
		DiscrimSignal = sum(respA .* RFWeight(:)');
		if (DiscrimSignal  > 0)
			LinPCorrect(cond) = LinPCorrect(cond) + 1;
		end
		if (DiscrimSignal == 0)
			LinPCorrect(cond) = LinPCorrect(cond) + 0.5;
		end
		DiscrimSignal = sum(respB) .* RFWeight(:);
		if (DiscrimSignal  < 0)
			LinPCorrect(cond) = LinPCorrect(cond) + 1;
		end
		if (DiscrimSignal == 0)
			LinPCorrect(cond) = LinPCorrect(cond) + 0.5;
		end
	end
	LinPCorrect(cond) = LinPCorrect(cond) / (2*NumTrials);
	NonLinPCorrect(cond) = NonLinPCorrect(cond) / (2*NumTrials);
	fprintf(1, 'Condition %d Probability Correct %d %d\n', cond, LinPCorrect(cond), NonLinPCorrect(cond));
end

DiscrimParams.MeanPhotoisoms = MeanPhotoisoms;
DiscrimParams.LinPCorrect = LinPCorrect;
DiscrimParams.NonlinPCorrect = NonLinPCorrect;
