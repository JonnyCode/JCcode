function [ReturnedCondition] = AverageAndVariance(CellInfo, Condition)
% AverageAndVariance.m
%	[ReturnedCondition] = AverageAndVariance(CellInfo, Condition)
%
%	Compute average and variance for all epochs in particular EpochCondition
%
%	Created: FMR sometime in the dark ages
%	Revised: 9/6/01 FMR and GDF
%		- updated to handle new CellInfo structures
%	    9/11/01 GDF
%		- updated to compute average and variance on matrices
%       12/08 FMR
%       - included standard usage and test
%
%%%SU 
%   load AverageAndVarianceTest/test
%   NewCondition = AverageAndVariance(CellInfo, OriginalCondition)
%
%%%TS isequal(NewCondition.AverageResponse, mean(OriginalCondition.EpochData.Data))
%%%TS isequal(NewCondition.VarianceResponse, var(OriginalCondition.EpochData.Data))

ReturnedCondition = Condition;
FieldFlag = isfield(Condition, 'EpochData');
ReturnedCellInfo = CellInfo;
ReturnedCondition = Condition;
NumEpochs = length(Condition.ExcludeEpochs);	
EpochPts = FindEpochPts(CellInfo, Condition);
[EpochData, Offset] = GetGoodEpochData(CellInfo, Condition);

if (size(EpochData, 1) > 1)
    ReturnedCondition.AverageResponse = mean(EpochData);
    ReturnedCondition.VarianceResponse = var(EpochData);
else
    ReturnedCondition.AverageResponse = EpochData;
    ReturnedCondition.VarianceResponse = nan;
end