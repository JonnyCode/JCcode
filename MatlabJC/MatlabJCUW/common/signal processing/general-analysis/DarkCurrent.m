function DarkCurrents = DarkCurrent(CellInfo, Condition)
%
%	DarkCurrents = DarkCurrent(CellInfo, Condition)
%
% This function estimates the dark current through time from the 
% saturating condition.
% Created: FMR ?
% Modified: GDR: 11/06/01
%			Now it handles new data format and works on matrices

NumEpochs = length(Condition.ExcludeEpochs);	
[EpochData, Offset] = GetGoodEpochData(CellInfo, Condition);
DarkCurrents = min(EpochData,[],2) + Offset';
