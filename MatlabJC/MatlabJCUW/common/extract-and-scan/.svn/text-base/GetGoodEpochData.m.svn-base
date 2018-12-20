function [GoodEpochData, Offset] = GetGoodEpochData(CellInfo, Condition);
%  GetGoodEpochData.m
%
%	[GoodEpochData, Offset] = GetGoodEpochData(CellInfo, Condition);
%
%  This funciton takes a CellInfo, and a Condition.  It then extracts the data
%  for the condition into a matrix.  It checks this against the ExcludeEpochs
%  vector and generates a matrix with only the non excluded epochs.

NumEpochs = length(Condition.ExcludeEpochs);
EpochFlag = 0;
[EpochData, TempOffset] = GetEpochData(CellInfo, Condition);
for epoch = 1:NumEpochs
	if (Condition.ExcludeEpochs(epoch) == 0); 
		EpochFlag = EpochFlag + 1;
		GoodEpochsList(EpochFlag) = epoch;
	end
end
GoodEpochData = EpochData(GoodEpochsList,:);
Offset = TempOffset(GoodEpochsList);
