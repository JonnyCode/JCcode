function Discriminant = FindDiscriminant2(CellInfo, RespConditionA, RespConditionB, TestCondition, FreqCutoff, RespLength, EpochList)
%
% FindDisciminant.m
%
%	Discriminant = FindDiscriminant(CellInfo, RespConditionA, RespConditionB, TestCondition, FreqCutoff, RespLength, EpochList)
%
%	Find mean discriminant for two alternative forced choice analysis.  Discriminant
%	returned is difference in mean of A and B responses.  All data (i.e. across all flash strengths)
%	except the "test" data set is used in discriminant calculation.
%
%	Created 10/00 FMR
%
%	Revisions
%		10/22/00 FMR removed mean correction
%		10/08/01 GDF made it compatible with new data format

% initialize space for discriminant
EpochPts = length(RespConditionA(TestCondition).AverageResponse);
Discriminant(1:EpochPts) = 0;
FlashPtA = round(RespConditionA(TestCondition).UserInfo.StimulusStartTime(1) / (FindSearchPara(RespConditionA(TestCondition), 'SampInterv') * 1e-6 * RespConditionA(TestCondition).DecimatePts));
FlashPtB = round(RespConditionB(TestCondition).UserInfo.StimulusStartTime(1) / (FindSearchPara(RespConditionB(TestCondition), 'SampInterv') * 1e-6 * RespConditionB(TestCondition).DecimatePts));

% loop across each flash strength
for cond = 1:length(RespConditionA)
	% if this is test condition only use data in EpochList for discriminant
	if (cond ~= TestCondition)
		[EpochDataA, Offset] = GetGoodEpochData(CellInfo, RespConditionA(cond));
		[EpochDataB, Offset] = GetGoodEpochData(CellInfo, RespConditionB(cond));
		AveRespA = sum(EpochDataA);
		AveRespB = sum(EpochDataB);
		Discriminant = Discriminant + AveRespA - AveRespB;
	end
end

% normalize to reasonable peak amplitude
Discriminant = Discriminant ./ abs(max(Discriminant) - min(Discriminant));

% set all frequencies above cutoff to zero
EpochPtsA = FindEpochPts(CellInfo, RespConditionA(TestCondition));
if (FreqCutoff > 0)
	Discriminant = fft(Discriminant);
	RawSamplingInterval = FindSearchPara(RespConditionA(TestCondition), 'SampInterv') * 1e-6;
	SamplingInterval = RawSamplingInterval * RespConditionA(TestCondition).DecimatePts;
	DeltaFreq = 1/(SamplingInterval * EpochPtsA);
	StartCutoff = floor(FreqCutoff/DeltaFreq);
	Discriminant(StartCutoff:EpochPtsA - StartCutoff) = 0;
	Discriminant = real(ifft(Discriminant));
end
Discriminant = Discriminant / (-min(Discriminant));

clear DiscriminantFilter;
DiscriminantFilter(1:length(Discriminant)) = 0;

