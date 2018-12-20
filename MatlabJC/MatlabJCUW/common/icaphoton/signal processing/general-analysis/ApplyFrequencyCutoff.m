function [ReturnedCellInfo, ReturnedCondition] = ApplyFrequencyCutoff(CellInfo, Condition, FreqCutoff);
% ApplyFrequencyCutoff.m
%
%  [ReturnedCellInfo, ReturnedCondition] = ApplyFrequencyCutoff(CellInfo, Condition, FreqCutoff)
%
% This function needs a CellInfo, a Condition and a FreqCutoff.  It then applies the 
% frequency cutoff to the data in CellInfo unless the Condition has a EpochData field
% inside the structure.  It returns a CellInfo and a Condition, with the data changed in
% the appropriate location.
%  Created: FMR  when?
%  Revised: GDF 09/06/01
%		    FMR 5/26/02
%				changed definition of sampling interval 

%  Get the Sampling interval from the search parameters for the condition 
clear ReturnedCondition ReturnedCellInfo FFTData EpochData
ReturnedCondition = Condition;
ReturnedCellInfo = CellInfo;
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
FFTData(:,FreqCutoffPts:length(FFTData(1,:))-FreqCutoffPts) = 0;
EpochData = real(ifft(FFTData, [], 2));
if (FieldFlag == 0)
	for epoch = 1:NumEpochs
		ReturnedCellInfo.EpochData.Data{Condition.EpochNumbers(epoch) + 1} = EpochData(epoch,:);
	end
end
if (FieldFlag == 1)
	ReturnedCondition.EpochData.Data = EpochData;
end

	
