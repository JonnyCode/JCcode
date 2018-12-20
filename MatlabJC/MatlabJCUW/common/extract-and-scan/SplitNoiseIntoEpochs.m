function ReturnedCondition = SplitNoiseIntoEpochs(CellInfo, Condition, EpochPts)
%
% ReturnCondition = SplitNoiseIntoEpochs(CellInfo, Condition, NoistPts)
%
%  Takes dark noise data from CellInfo and breaks it up into the Epochs
%  that are EpochPts in length
%  Created  FMR   ?
%  Revised GDF  9/14/01
%			- it now gets the data from CellInfo and is compatible with new data format

NumEpochs = length(Condition.ExcludeEpochs) - sum(Condition.ExcludeEpochs);		% Get the number of epochs
NoisePts = FindEpochPts(CellInfo, Condition);
Removed = {['EpochNumbers'], ['AverageResponse'], ['VarianceResponse'], ['EpochData']};
ReturnedCondition = rmfield(Condition, Removed);
ReturnedCondition.Label = 'Noise_Epochs';
ReturnedCondition.UserInfo.EpochPts = EpochPts;

NumNewEpochs = floor(NoisePts / EpochPts);
ReturnCondition.UserInfo.StimulusStartTime(1:NumNewEpochs*length(Condition.ExcludeEpochs)) = 0;

clear ReturnedCondition.ExcludeEpochs;
[EpochData, Offset] = GetGoodEpochData(CellInfo, Condition);

% go through each epoch and combine pairs of vectors
for cnt = 1:NumEpochs
	for epoch = 1:NumNewEpochs
		ReturnedCondition.EpochData.Data((cnt-1) * NumNewEpochs + epoch, :) = EpochData(cnt, (epoch-1)*EpochPts+1: epoch*EpochPts);
		ReturnedCondition.EpochData.Offset((cnt-1) * NumNewEpochs + epoch) = Offset(cnt);
		ReturnedCondition.ExcludeEpochs((cnt-1) * NumNewEpochs + epoch) = 0;
	end
end
