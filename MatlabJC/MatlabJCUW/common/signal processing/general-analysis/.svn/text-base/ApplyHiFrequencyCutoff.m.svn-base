function ReturnedCondition = ApplyHiFrequencyCutoff(CellInfo, Condition, FreqCutoff);
% ApplyHiFrequencyCutoff.m
%
%  [ReturnedCellInfo, ReturnedCondition] = ApplyFrequencyCutoff(CellInfo, Condition, FreqCutoff)
%
% Hi pass filter condition.  Routine sets all frequency components below specified cutoff to 0.
%
%	11/03 Created FMR
%

% Get the Sampling interval from the search parameters for the condition 
clear ReturnedCondition FFTData EpochData
ReturnedCondition = Condition;
RawSamplingInterval = FindSearchPara(Condition, 'SampInterv');
SamplingInterval = RawSamplingInterval * Condition.DecimatePts * 1e-6;
EpochPts = FindEpochPts(CellInfo, Condition);
FieldFlag = isfield(Condition, 'EpochData');

FreqStepSize = 1/(SamplingInterval * EpochPts);
FreqCutoffPts = round(FreqCutoff / FreqStepSize);
NumEpochs = length(Condition.ExcludeEpochs);

% Apply the frequency cutoff to each epoch specified in the EpochNumbers of EpochCondition
[EpochData, Offset] = GetEpochData(CellInfo, Condition);
% set upper frequency limit and decimate
FFTData = fft(EpochData, [], 2);
if (FreqCutoffPts > 1)
	for epoch = 1:size(FFTData, 1)
		FFTData(epoch,1:FreqCutoffPts) = 0;
		FFTData(epoch,EpochPts-FreqCutoffPts:EpochPts) = 0;
	end
end
EpochData = real(ifft(FFTData, [], 2));
ReturnedCondition.EpochData.Data = EpochData;
ReturnedCondition = AverageAndVariance(CellInfo, ReturnedCondition);
	
