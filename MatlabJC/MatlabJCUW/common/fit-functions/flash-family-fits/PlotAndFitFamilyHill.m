function [CellParameters, PlotData] = PlotAndFitFamilyHill(FamilyCondition, FamCond, Polarity)

clf; subplot(1, 2, 1); hold on
clear FlashStrength;
NumFlash = 1;

if (isfield(FamilyCondition(FamCond).UserInfo, 'FamilyStep'))
    FamilyStep = FamilyCondition(FamCond).UserInfo.FamilyStep;
else
    FamilyStep = FamilyCondition(FamCond).FamilyStep;
    fprintf(1, 'Warning: no calibration - using uncalibrated amps\n');
end
FlashStrength = unique(FamilyStep);
FamilyStep = FamilyCondition(FamCond).FamilyStep;
UniqueFamilyStep = unique(FamilyStep);
tme = 1:length(FamilyCondition(FamCond).EpochData.Data);
tme = tme * FindSearchPara(FamilyCondition(FamCond), 'SampInterv') * FamilyCondition(FamCond).DecimatePts * 1e-6;
for cond = 1:length(UniqueFamilyStep)
    [EpochData, Offset] = ExtractFamilyData(FamilyCondition(FamCond), UniqueFamilyStep(cond));
    if (size(EpochData, 1) > 1)
        PlotData.Data(cond, :) = mean(EpochData);
    else
        PlotData.Data(cond, :) = EpochData;
    end
    plot(tme, PlotData.Data(cond, :));
    if (Polarity == -1)
        response(NumFlash) = -min(PlotData.Data(cond, :));
    else
        response(NumFlash) = max(PlotData.Data(cond, :));
    end
    NumFlash = NumFlash + 1;
end

xlabel('time');
ylabel('response');
hold off;
subplot(1, 2, 2);
maxresponse = max(response);
response = response / maxresponse;

plot(log10(FlashStrength), response, 'o');
xlabel('log(Flash Strength)');
ylabel('R/R_{max}');

% fit with model
coef = {1 1};
coef = cat(1, coef{:});
fitcoef = nlinfit(FlashStrength, log10(response), 'hill', coef);
bestfit = hill(fitcoef, FlashStrength);
plot(0);
hold  on
plot(log10(FlashStrength), log10(response), 'o');
plot(log10(FlashStrength), bestfit);
xlabel('log(Flash Strength)');
ylabel('log(R/R_{max})');
hold off
fprintf('\toffset = %d\tpower = %d\n', fitcoef(1), fitcoef(2));
CellParameters.HalfSatIntensity = fitcoef(1);
CellParameters.MaxResponse = maxresponse;

PlotData.response = response;
PlotData.FlashStrength = FlashStrength;
PlotData.Fit = bestfit;
