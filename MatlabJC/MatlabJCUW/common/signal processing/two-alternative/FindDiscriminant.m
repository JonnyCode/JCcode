function Discriminant = FindDiscriminant(CellInfo, RespConditionA, RespConditionB, TestCondition, FreqCutoff, RespLength, EpochList)
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
	if (cond == TestCondition)
		[EpochDataA, OffsetA] = GetGoodEpochData(CellInfo, RespConditionA(cond));
		[EpochDataB, OffsetB] = GetGoodEpochData(CellInfo, RespConditionB(cond));
		% make sure valid epoch
		GoodEpochs = 0;
		for cnt = 1:length(EpochList)
			if (RespConditionA(cond).ExcludeEpochs(EpochList(cnt)) == 0)
				if (RespConditionB(cond).ExcludeEpochs(EpochList(cnt)) == 0)
					GoodEpochs = GoodEpochs + 1;
					TempEpochList(GoodEpochs) = EpochList(cnt);
				end
			end
		end
		if (length(TempEpochList) > 0)
			clear EpochList
			EpochList = TempEpochList;			
			AveRespA = sum(EpochDataA(EpochList,:), 1);				
			AveRespB = sum(EpochDataB(EpochList,:), 1);
			Discriminant = Discriminant + AveRespA - AveRespB;
		end
		% if not test condition use all data
	else	
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

% truncate discriminant if needed
%if (RespLength < EpochPtsA)
if (0)
	fprintf(1, 'truncating discriminant\n');
	FlashA = RespConditionA(TestCondition).UserInfo.StimulusStartTime(1) / SamplingInterval;
	if ((FlashA + RespLength) > EpochPtsA)
		EndPoint = EpochPtsA;
	else
		EndPoint = FlashA + RespLength;
	end
	DiscriminantFilter(FlashA+1:EndPoint) = 1;
	EpochPtsB = FindEpochPts(CellInfo, RespConditionB(TestCondition));
	FlashB = RespConditionB(TestCondition).UserInfo.StimulusStartTime(1) / SamplingInterval;
	if ((FlashB + RespLength) > EpochPtsB)
		EndPoint = EpochPtsB;
	else
		EndPoint = FlashB + RespLength;
	end
	DiscriminantFilter(FlashB+1:EndPoint) = 1;
	Discriminant = Discriminant .* DiscriminantFilter;
end
