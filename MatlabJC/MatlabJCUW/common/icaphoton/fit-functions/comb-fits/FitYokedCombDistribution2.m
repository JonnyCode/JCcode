function [CombParameters, CombFit] = FitYokedCombDistribution2(RespProjection)

%	FitYokedCombDistribution2
%
%	Fit comb distribution (i.e. sum of gaussians with weights determined from Poisson dist)
%	to a series of histograms with flash strengths specified in MeanPhotoisoms.  RespProjection
% 	should contain a noise histogram (NoiseHist) and histograms for A and B responses (RespAHist
%	and RespBHist).  
%

global DarkNoiseSD;
global FlashStrength;
global NumResponses;
global NumHistBins;
global NumberCombConditions;

clear Respx;
Respx(1:NumHistBins) = -NumHistBins/2:NumHistBins/2-1;
Respx = (floor(max(RespProjection(NumberCombConditions).RespAProjection)) + 1) * 2 * Respx / NumHistBins;
Respx = Respx';

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

coef = {0.3 1 1};
coef = cat(1, coef{:});
	
fitcoef = nlinfit(YokedRespx, YokedRespHist, 'yoked_comb_fit2', coef);
CombFit.FitA = yoked_comb_fit2(fitcoef, YokedRespx);

CombParameters.DarkNoiseSD = DarkNoiseSD;
CombParameters.UnitaryAmpA = fitcoef(2);
CombParameters.UnitarySDA = fitcoef(1);
CombParameters.CollectingArea = fitcoef(3);

