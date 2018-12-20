function [NewCondition, BadEpochs] = ScanEpochs(Condition, BadEpochs)
%
%	function [NewCondition, BadEpochs] = ScanEpochs(Condition, BadEpochs)
%
%	This function plots each epoch individually and allows the user
%	to mark certain epochs for exclusion.

[NumEpochs, EpochPts] = size(Condition.EpochData.Data);	
NewCondition = Condition;
tme = 1:EpochPts;
epoch = 1;

while (epoch <= NumEpochs)
	if (Condition.ExcludeEpochs(epoch) == 1)
		NewCondition.ExcludeEpochs(epoch) = 1;
	end
	if (Condition.ExcludeEpochs(epoch) == 0)
		tme = 1:length(Condition.EpochData.Data(epoch, :));
		plot(tme, Condition.EpochData.Data(epoch, :), tme, Condition.AverageResponse);
%		plot(tme, Condition.EpochData(epoch, :));
		cnt = 0;
		ExcludeFlag = input('Exclude?');
		if (ExcludeFlag == 1)
			NewCondition.ExcludeEpochs(epoch) = 1;
			BadEpochs(length(BadEpochs) + 1) = Condition.EpochNumbers(epoch);
		end
		if (ExcludeFlag == 2)
			break;
		end
		if (ExcludeFlag == -1)
			epoch = epoch - 2;
		end
	end
	fprintf(1, 'epoch #%d\n', epoch);
	epoch = epoch + 1;
end
