function ReturnCondition = BaselineCorrect(CellInfo, Condition, StartBaseline, EndBaseline)
%
%	function ReturnCondition = BaselineCorrect(CellInfo, Condition, StartBaseline, EndBaseline)
%
%	This function takes a small junk of baseline and subtracts off the average of that
%	baseline from the rest of the data.  This must always be applied to a condition structure
%
%	Created: GDF 11/01

[NumEpochs, EpochPts] = size(Condition.EpochData.Data);
ReturnCondition = Condition;

for epoch = 1:NumEpochs
	temp(1:EpochPts) = Condition.EpochData.Data(epoch, 1:EpochPts);
	baseline(1:EndBaseline-StartBaseline) = Condition.EpochData.Data(epoch, (StartBaseline+1):EndBaseline);
	temp = temp - (sum(baseline) / (EndBaseline - StartBaseline));
	ReturnCondition.EpochData.Data(epoch, 1:EpochPts) = temp(:);
end
