function FamilyCondition = LoadAndSmoothFamilyCondition(CellInfo, FreqCutoff, DecimatePts)

clear FamilyCondition;
FamilyCondition = CellInfo.FamilyCondition;
for cond = 1:length(FamilyCondition)
	[EpochData, Offset] = GetEpochData(CellInfo, CellInfo.FamilyCondition(cond));
	FamilyCondition(cond).EpochData.Data = EpochData;
	FamilyCondition(cond).EpochData.Offset = Offset;
	if (DecimatePts > 1)
		FamilyCondition(cond) = DecimateCondition(CellInfo, FamilyCondition(cond), DecimatePts);
	end
end

% smoothing: FamilyCondition
for cond = 1:length(FamilyCondition)
	[CellInfo, Condition] = ApplyFrequencyCutoff(CellInfo, FamilyCondition(cond), FreqCutoff);
	FamilyCondition(cond) = Condition;
end
