function [CombParameters, CombFit] = FitYokedCombDistSingles(RespProjection)
%
%	function [CombParameters, CombFit] = FitYokedCombDistSingles(RespProjection)
%
%	Created: FMR ?
%	Revised: GDF made to extract singles.


global DarkNoiseSD;
global MeanPhotoisoms;
global NumResponses;
global NumberCombConditions;
global NumHistBins;

clear Respx;
Respx = RespProjection(1).Respx;

% combine and fit noise distributions
NoiseHist(1:NumHistBins) = 0;
NumResponses = 0;
for cond = 1:NumberCombConditions
	NoiseHist = NoiseHist + RespProjection(cond).NoiseHist;
	NumResponses = NumResponses + length(RespProjection(cond).NoiseProjection);
end

coef = {1};
coef = cat(1, coef{:});
	
fitcoef = nlinfit(Respx, NoiseHist, 'ZeroMeanGaussian', coef);
fit = ZeroMeanGaussian(fitcoef, Respx);
plot(Respx, NoiseHist, Respx, fit);
DarkNoiseSD = fitcoef(1);

% fit A responses
clear YokedRespHist;
clear YokedRespx;
YokedRespHist(1:NumberCombConditions*NumHistBins) = 0;
YokedRespx(1:NumberCombConditions*NumHistBins) = 0;
YokedRespx = YokedRespx';
YokedRespHist = YokedRespHist';

for cond = 1:NumberCombConditions
	YokedRespHist((cond-1)*NumHistBins+1:cond*NumHistBins) = RespProjection(cond).RespAHist;	
	YokedRespx((cond-1)*NumHistBins+1:cond*NumHistBins) = Respx(1:NumHistBins);	
	NumResponses(cond) = length(RespProjection(cond).RespAProjection);
end

coef = {0.3 1};
coef = cat(1, coef{:});
	
fitcoef = nlinfit(YokedRespx, YokedRespHist, 'yoked_comb_fit', coef);
CombFit.FitA = yoked_comb_fit(fitcoef, YokedRespx);

CombParameters.DarkNoiseSD = DarkNoiseSD;
CombParameters.UnitaryAmpA = fitcoef(2);
CombParameters.UnitarySDA = fitcoef(1);

% fit B responses
%clear YokedRespHist;
%YokedRespHist(1:NumberCombConditions*NumHistBins) = 0;
%YokedRespHist = YokedRespHist';

%for cond = 1:NumberCombConditions
%	YokedRespHist((cond-1)*NumHistBins+1:cond*NumHistBins) = RespProjection(cond).RespBHist;	
%	YokedRespx((cond-1)*NumHistBins+1:cond*NumHistBins) = Respx(1:NumHistBins);	
%	NumResponses(cond) = length(RespProjection(cond).RespBProjection);
%end
%
%coef = {0.3 1};
%coef = cat(1, coef{:});
%	
%fitcoef = nlinfit(YokedRespx, YokedRespHist, 'yoked_comb_fit', coef);
%CombFit.FitB = yoked_comb_fit(fitcoef, YokedRespx);
%
%CombParameters.DarkNoiseSD = DarkNoiseSD;
%CombParameters.UnitaryAmpB = fitcoef(2);
%CombParameters.UnitarySDB = fitcoef(1);


