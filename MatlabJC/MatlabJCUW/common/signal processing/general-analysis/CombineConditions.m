function ReturnCondition = CombineConditions(CellInfo, ConditionA, ConditionB)
% CombineConditions
%
% Combine two conditions into one.
%
%   Revised 7/18/05 FMR - checked whether EpochData defined and if not
%   combined without data

NumEpochs = length(ConditionB.ExcludeEpochs);	
ReturnCondition = ConditionA;
StartEpoch = length(ConditionA.ExcludeEpochs) + 1;

ReturnCondition.EpochNumbers(StartEpoch:StartEpoch + NumEpochs - 1) = ConditionB.EpochNumbers(:);
ReturnCondition.ExcludeEpochs(StartEpoch:StartEpoch + NumEpochs - 1) = ConditionB.ExcludeEpochs(:);
ReturnCondition.EpochTime(StartEpoch:StartEpoch + NumEpochs - 1) = ConditionB.EpochTime(:);

if (isfield(ConditionA.UserInfo, 'StimulusStartTime'))
    ReturnCondition.UserInfo.StimulusStartTime(StartEpoch:StartEpoch + NumEpochs - 1) = ConditionB.UserInfo.StimulusStartTime(:);
end

if (isfield(ConditionA, 'EpochData'))
    for epoch = 1:NumEpochs
        ReturnCondition.EpochData.Data(StartEpoch + epoch - 1, :) = ConditionB.EpochData.Data(epoch, :);
        ReturnCondition.EpochData.Offset(StartEpoch + epoch - 1) = ConditionB.EpochData.Offset(epoch);
    end
    ReturnCondition = AverageAndVariance(CellInfo, ReturnCondition);
end

