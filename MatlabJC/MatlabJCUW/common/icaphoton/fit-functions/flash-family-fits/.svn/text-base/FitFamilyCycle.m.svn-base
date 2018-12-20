function [HalfSatFlash, FamilyTime] = FitFamilyCycle(CellInfo, FamilyCondition)

clear FlashStrength;
if (isfield(FamilyCondition.UserInfo, 'FamilyStep'))
	FamilyStep = FamilyCondition.UserInfo.FamilyStep;
else
	FamilyStep = FamilyCondition.FamilyStep;
	fprintf(1, 'Warning: no calibration - using uncalibrated amps\n');
end
FlashStrength = unique(FamilyStep)
FamilyStep = FamilyCondition.FamilyStep;
UniqueFamilyStep = unique(FamilyStep);
NumFamilyRepeats = length(FamilyStep) / length(UniqueFamilyStep);
[EpochData, Offset] = GetEpochData(CellInfo, FamilyCondition);
ResponseAmps = -min(EpochData');
coef = [0.001];
for FamCond = 1:NumFamilyRepeats-1
	for cond = 1:length(UniqueFamilyStep)
		response(cond) = ResponseAmps(FamCond*length(UniqueFamilyStep) + cond);
	end
	maxresponse = max(response);
	response = response / maxresponse;

	% fit with model
	fitcoef = nlinfit(FlashStrength', response', 'familymodel', coef);
    fit = familymodel(fitcoef, FlashStrength');
	HalfSatFlash(FamCond) =(1/fitcoef(1))*log(2);
    FamilyTime(FamCond) = FamilyCondition.EpochTime(FamCond * length(UniqueFamilyStep));
end

