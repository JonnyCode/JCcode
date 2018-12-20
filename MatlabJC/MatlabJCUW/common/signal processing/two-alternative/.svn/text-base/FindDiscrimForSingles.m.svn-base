function Discriminant = FindDiscrimForSingles(CellInfo, RespConditionA, NoiseCondition, TestCondition, FreqCutoff, RespLength, EpochList)
%
% FindDisciminant.m
%
%	Discriminant = FindDiscriminantForSingles(CellInfo, RespConditionA, NoiseCondition, TestCondition, FreqCutoff, RespLength, EpochList)
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
%		12/01	GDF revised to extract singles

% initialize space for discriminant
EpochPts = FindEpochPts(CellInfo, RespConditionA(TestCondition));
Discriminant(1:EpochPts) = 0;

% loop across each flash strength
for cond = 1:length(RespConditionA)
	clear TempEpochList
	% if this is test condition only use data in EpochList for discriminant
	if (cond == TestCondition)
		[EpochDataA, OffsetA] = GetGoodEpochData(CellInfo, RespConditionA(cond));
		GoodEpochs = 0;
		for epoch = 1:(length(NoiseCondition.ExcludeEpochs))
			if (NoiseCondition.ExcludeEpochs(epoch) == 0)
				GoodEpochs = GoodEpochs + 1;
				EpochDataNoise(GoodEpochs,:) = NoiseCondition.EpochData.Data(epoch,:);
			end
		end
		% make sure valid epoch
		GoodEpochs = 0;
		for cnt = 1:length(EpochList)
			if (RespConditionA(cond).ExcludeEpochs(EpochList(cnt)) == 0)
				if (NoiseCondition.ExcludeEpochs(EpochList(cnt)) == 0)
					GoodEpochs = GoodEpochs + 1;
					TempEpochList(GoodEpochs) = EpochList(cnt);
				end
			end
		end
		clear EpochList
		EpochList = TempEpochList;			
		AveRespA = mean(EpochDataA(EpochList,:), 1);				
		AveRespB = mean(EpochDataNoise(EpochList,:), 1);
		Discriminant = Discriminant + AveRespA - AveRespB;
		clear AveRespA AveRespB
	% if not test condition use all data
	else	
		[EpochDataA, Offset] = GetGoodEpochData(CellInfo, RespConditionA(cond));
		GoodEpochs = 0;
		for epoch = 1:(length(NoiseCondition.ExcludeEpochs))
			if (NoiseCondition.ExcludeEpochs(epoch) == 0)
				GoodEpochs = GoodEpochs + 1;
				EpochDataNoise(GoodEpochs,:) = NoiseCondition.EpochData.Data(epoch,:);
			end
		end		% make sure valid epoch
		AveRespA = mean(EpochDataA, 1);
		AveRespNoise = mean(EpochDataNoise, 1);
		Discriminant = Discriminant + AveRespA - AveRespNoise;
	end
end

% normalize to reasonable peak amplitude
Discriminant = Discriminant ./ abs(max(Discriminant) - min(Discriminant));

% set all frequencies above cutoff to zero
if (FreqCutoff > 0)
	Discriminant = fft(Discriminant);
	RawSamplingInterval = 1 / FindSearchPara(RespConditionA(TestCondition), 'SampInterv');
	SamplingInterval = RawSamplingInterval * RespConditionA(TestCondition).DecimatePts;
	EpochPtsA = FindEpochPts(CellInfo, RespConditionA(TestCondition));
	DeltaFreq = 1/(SamplingInterval * EpochPtsA);
	StartCutoff = floor(FreqCutoff/DeltaFreq);
	Discriminant(StartCutoff:EpochPtsA - StartCutoff) = 0;
	Discriminant = real(ifft(Discriminant));
end
Discriminant = Discriminant / (-min(Discriminant));


clear DiscriminantFilter;
DiscriminantFilter(1:length(Discriminant)) = 0;

% truncate discriminant if needed
if (RespLength < EpochPtsA)
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

plot(Discriminant);
