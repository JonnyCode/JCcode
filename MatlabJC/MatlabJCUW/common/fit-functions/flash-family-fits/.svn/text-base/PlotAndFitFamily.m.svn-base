function [CellParameters, PlotData] = PlotAndFitFamily(FamilyCondition, EpochCondition, CellType, MinFlag)

NumFamCond = length(FamilyCondition);
subplot(1, 2, 1)
plot(0);
hold on;
clear FlashStrength;
NumFlash = 1;
for FamCond = 1:1 %NumFamCond
	if (CellType == 'rod')
        if (isfield(FamilyCondition(FamCond).UserInfo, 'RodStimulusStrength'))
            FamilyStep = FamilyCondition(FamCond).UserInfo.RodStimulusStrength;
        else
            FamilyStep = FamilyCondition(FamCond).FamilyStep;
            fprintf(1, 'Warning: no calibration - using uncalibrated amps\n');
        end
    else
        if (isfield(FamilyCondition(FamCond).UserInfo, 'ConeStimulusStrength'))
            FamilyStep = FamilyCondition(FamCond).UserInfo.ConeStimulusStrength;
        else
            FamilyStep = FamilyCondition(FamCond).FamilyStep;
            fprintf(1, 'Warning: no calibration - using uncalibrated amps\n');
        end
    end
    
    if (FamCond == 1)
		FlashStrength = unique(FamilyStep);
		RespLength = size(FamilyCondition(FamCond).EpochData.Data, 2);
	else
		FlashStrength = [FlashStrength, unique(FamilyStep)];
	end
	FamilyStep = unique(FamilyCondition(FamCond).FamilyStep);
	tme = 1:length(FamilyCondition(FamCond).EpochData.Data);
	tme = tme * FindSearchPara(FamilyCondition(FamCond), 'SampInterv') * FamilyCondition(FamCond).DecimatePts * 1e-6;
	for cond = 1:length(unique(FamilyStep))
		EpochData = ExtractFamilyData(FamilyCondition(FamCond), FamilyStep(cond));
        if (cond == 1)
            NumPts = size(EpochData, 2);
        end
        if (MinFlag)
            PlotData.Data(NumFlash, 1:NumPts) = -mean(EpochData(:, 1:NumPts));
            response(NumFlash) = max(-mean(EpochData));
    		plot(tme, -mean(EpochData));
        else
            PlotData.Data(NumFlash, 1:NumPts) = mean(EpochData(:, 1:NumPts));
            response(NumFlash) = max(mean(EpochData));
    		plot(tme, mean(EpochData));
        end
    	NumFlash = NumFlash + 1;
	end
end
if (isstruct(EpochCondition))
    NumCond = length(EpochCondition);
    for cond = 1:NumCond
        if (strncmp(EpochCondition(cond).Label, 'dim', 3))
            NumPts2 = length(EpochCondition(cond).AverageResponse);
            tme = 1:NumPts2;
            tme = tme * FindSearchPara(EpochCondition(cond), 'SampInterv') * EpochCondition(cond).DecimatePts * 1e-6;
            if (MinFlag)
                PlotData.Data(NumFlash, :) = -EpochCondition(cond).AverageResponse(1:NumPts);
                response(NumFlash) = max(-EpochCondition(cond).AverageResponse);
                plot(tme, -EpochCondition(cond).AverageResponse);
            else
                PlotData.Data(NumFlash, 1:NumPts) = EpochCondition(cond).AverageResponse(1:NumPts);
                response(NumFlash) = max(EpochCondition(cond).AverageResponse);
                plot(tme, EpochCondition(cond).AverageResponse);
            end
            NumFlash = NumFlash + 1;
            FlashStrength = [FlashStrength, EpochCondition(cond).UserInfo.StimulusAmp];
        end
    end
end

xlabel('time');
ylabel('response');
axis tight;
hold off;
subplot(1, 2, 2);
maxdarkcurrent = max(response);
response = response / maxdarkcurrent;

plot(log10(FlashStrength), response, 'o');
xlabel('log(Flash Strength)');
ylabel('R/R_{max}');


% fit with model
coef = {0.001};
coef = cat(1, coef{:});
fitcoef = nlinfit(FlashStrength', response', 'familymodel', coef);
fit = familymodel(fitcoef, FlashStrength);
plot(0);
hold  on
plot(log10(FlashStrength), response, 'o');
plot(log10(FlashStrength), fit);
HalfSatIntensity =(1/fitcoef(1))*log(2);
CellParameters.HalfSatIntensity = HalfSatIntensity;
axis tight;
hold off

PlotData.response = response;
PlotData.FlashStrength = FlashStrength;
PlotData.Fit = fit;
PlotData.NumCycles = length(FamilyStep) / length(unique(FamilyStep));
PlotData.DarkCurrent = maxdarkcurrent;
