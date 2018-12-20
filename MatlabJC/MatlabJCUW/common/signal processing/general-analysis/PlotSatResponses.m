function PlotSatResponses(FamilyCondition, label)
% PlotSatResponse
% 
% Plot amplitude of saturating responses from flash families to check stability.

NumFamCond = length(FamilyCondition);
NumFlash = 1;
for FamCond = 1:NumFamCond
	if (strcmp(FamilyCondition(FamCond).Label, label) == 1)
		FamilyStep = FamilyCondition(FamCond).FamilyStep;
		Indices = find(FamilyStep == max(FamilyStep));
		EpochData = ExtractFamilyData(FamilyCondition(FamCond), max(FamilyStep));
		plot(Indices, min(EpochData') - max(EpochData'), 'o')
	end
end
