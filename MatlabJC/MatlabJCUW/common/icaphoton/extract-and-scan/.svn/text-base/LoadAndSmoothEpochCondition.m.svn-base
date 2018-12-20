function EpochCondition = LoadAndSmoothEpochCondition(CellInfo, FreqCutoff, DecimatePts);

clear EpochCondition;
EpochCondition = CellInfo.EpochCondition;
for cond = 1:length(EpochCondition)
	[EpochData, Offset] = GetEpochData(CellInfo, CellInfo.EpochCondition(cond));
	EpochCondition(cond).EpochData.Data = EpochData;
	EpochCondition(cond).EpochData.Offset = Offset;
	EpochCondition(cond) = AverageAndVariance(CellInfo, EpochCondition(cond));
	if (DecimatePts > 1)
		EpochCondition(cond) = DecimateCondition(CellInfo, EpochCondition(cond), DecimatePts);
	end
end

% smoothing: EpochCondition
for cond = 1:length(EpochCondition)
	[CellInfo, Condition] = ApplyFrequencyCutoff(CellInfo, EpochCondition(cond), FreqCutoff);
	Condition = AverageAndVariance(CellInfo, Condition);
	EpochCondition(cond) = Condition;
end
