function [fit] = yoked_comb_fit(beta, x)

MaxHit = 40;
global DarkNoiseSD;
global MeanPhotoisoms;
global NumResponses;
global NumHistBins;
global NumberCombConditions;

fit(1:length(x)) = 0;
fit = fit';
xinc = x(2) - x(1);
fitx(1:NumHistBins) = x(1:NumHistBins);
tempfit(1:NumHistBins) = 0;

for cond = 1:NumberCombConditions
	tempfit(1:NumHistBins) = 0;
	fitx(1:NumHistBins) = x((cond-1)*NumHistBins+1:cond*NumHistBins);		
	for n = 0:MaxHit
		% Poisson factor for probability of n photoisomerizations
		temp = (exp(-MeanPhotoisoms(cond)) * MeanPhotoisoms(cond)^n / factorial(n));
		% multiply by normalization factor for Gaussian
		temp = temp * ((2 * 3.141592 * (DarkNoiseSD^2 + n * beta(1)^2))^(-0.5));
		% multiply by Gaussian with mean n*single-amp and variance dark + n*single variance
		temp = temp .* exp(-((fitx - n * beta(2)).^2) / (2 * (DarkNoiseSD^2 + n * beta(1)^2)));
		% add to fit for all photon counts
		tempfit = tempfit + temp .* NumResponses(cond) .* xinc;
	end
	fit((cond-1)*NumHistBins+1:cond*NumHistBins) = tempfit(1:NumHistBins);
end
