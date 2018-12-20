function [EpochData, Offset] = ExtractFamilyData(FamilyCondition, FamilyStep)
%
%	function EpochData = ExtractFamilyData(Condition, BadEpochs)
%
%	This function extracts data from a FamilyCondition with a 
%	particular FamilyStep (e.g. amp) and returns it in a a matrix
%	EpochData with the epochs and a vector Offset with the offset
%	applied to each epoch when the data was read in.
%
%	Created ? FMR

[NumEpochs, EpochPts] = size(FamilyCondition.EpochData.Data);	
epoch = 1;
GoodEpochs = 1;

while (epoch <= NumEpochs)
	if (FamilyCondition.FamilyStep(epoch) == FamilyStep)
		if (FamilyCondition.ExcludeEpochs(epoch) == 0)
			EpochData(GoodEpochs, :) = FamilyCondition.EpochData.Data(epoch, :);
			Offset(GoodEpochs) = FamilyCondition.EpochData.Offset(epoch);
			GoodEpochs = GoodEpochs + 1;
		end
	end
	epoch = epoch + 1;
end
