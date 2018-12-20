function [ResponseA, ResponseB] = SampleFromCombDistribution(CombParameters, NumReceptors, NumTrials)
%
%	[ResponseA, ResponseB] = SampleFromCombDistribution(CombParameters, NumReceptors, NumTrials)
%
%	This is a revised version that is faster because it uses MatLabs
% 	ability to generate a matrix of samples rather than one at a time.
%	Created FMR ?
%	Revised GDF: Faster, samples one matrix at a time

if (CombParameters.RateA > 0)
	NumHits = poissrnd(CombParameters.RateA, NumTrials, NumReceptors);
else
	NumHits = 0;
end

% Calc standard devs
StDev = sqrt(NumHits * CombParameters.UnitarySDA.^2 + CombParameters.DarkNoiseSD.^2);

% Generate samples from a normal distribution and multiply by the stdev.
[trials, receptors] = size(NumHits);
NormDist = normrnd(0, 1, trials, receptors);
ResponseA = NormDist .* StDev;
ResponseA = ResponseA + (NumHits * CombParameters.UnitaryAmpA);

clear NumHits trials receptors

if (CombParameters.RateB > 0)
	NumHits = poissrnd(CombParameters.RateB, NumTrials, NumReceptors);
else
	NumHits = 0;
end

% Calc standard devs
StDev = sqrt(NumHits * CombParameters.UnitarySDB.^2 + CombParameters.DarkNoiseSD.^2);

% Generate samples from a normal distribution and multiply by the stdev.
[trials, receptors] = size(NumHits);
NormDist = normrnd(0, 1, trials, receptors);
ResponseB = NormDist .* StDev;
ResponseB = ResponseB + (NumHits * CombParameters.UnitaryAmpB);

clear NumHits trials receptors
