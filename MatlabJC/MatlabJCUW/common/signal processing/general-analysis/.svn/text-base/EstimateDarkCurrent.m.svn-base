function ReturnedCellInfo = EstimateDarkCurrent(CellInfo, SatConditionLabel)
%
% 	function ReturnedCondition = EstimateDarkCurrent(SatCondition, SatConditionLabel)
%
% Estimates the dark current during the experiment from the saturating flashes
% The SatConditionLabel should be a string that identifies the label you have 
% associated with your a Condition that delivers a saturating flash.
% Created FMR: ?
% Modified GDF: 11/06/01
%			This has been generalized from the last version.
%			It works with the new data format, and it finds the
%			Condition inside CellInfo containing the saturating flash and uses this
%			To compute the dark current over the course of the experiment.

ReturnedCellInfo = CellInfo;

% Find Condition with saturating flash
satcond = 0;
NumConds = length(CellInfo.EpochCondition);
for cond = 1:NumConds
	if (strcmp(CellInfo.EpochCondition(cond).Label, SatConditionLabel))
		satcond = cond;
	end
end
if satcond == 0
	fprintf(1, 'Warning: no saturating condition was found. \n');
end

% estimate dark current level over time from saturating responses
CellInfo.EpochCondition(satcond).UserInfo.DarkCurrent = DarkCurrent(CellInfo, CellInfo.EpochCondition(satcond));
SatCondition = CellInfo.EpochCondition(satcond);

% add dark current to offset slot in each epoch
for cond = 1:NumConds
	[EpochData, Offset] = GetEpochData(CellInfo, CellInfo.EpochCondition(cond));
	Condition = CellInfo.EpochCondition(cond);
	NumEpochs = length(Condition.ExcludeEpochs);	
	for epoch = 1:NumEpochs
		EpochTime = Condition.EpochTime(epoch);
		for satepoch = 1:length(SatCondition.UserInfo.DarkCurrent)
			if (SatCondition.EpochTime(satepoch) > EpochTime)
				if (SatCondition.ExcludeEpochs(satepoch) == 0)
					break
				else
					satepoch = satepoch - 1;
					break
				end
			end
		end
		WeightFact = (EpochTime - SatCondition.EpochTime(satepoch-1)) / (SatCondition.EpochTime(satepoch) - SatCondition.EpochTime(satepoch-1));
		DarkCurrent = SatCondition.UserInfo.DarkCurrent(satepoch-1) * (1 - WeightFact) + SatCondition.UserInfo.DarkCurrent(satepoch) * WeightFact;	
		Condition.UserInfo.DarkCurrent(epoch) = Offset(epoch) - DarkCurrent;
	end
	ReturnedCellInfo.EpochCondition(cond) = Condition;
end
