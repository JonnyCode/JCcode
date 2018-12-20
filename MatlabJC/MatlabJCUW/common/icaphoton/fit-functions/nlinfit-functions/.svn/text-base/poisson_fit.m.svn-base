% poisson_fit
%
% 	function to fit elementary response histogram with free parameters
% 	for baseline standard deviation, elementary response standard deviation,
%	mean number of events per trial, mean event amplitude, and total number of
%	trials.
%
%	Created 5/00 FMR

function [fit] = poisson_fit(beta, x)

fit(1:length(x)) = 0;
fit = fit';
xinc = x(2) - x(1);

MaxHit = 20;
global DarkNoiseSD;
global MeanPhotoisoms;
global NumResponses;

for n = 0:MaxHit
	% Poisson factor for probability of n photoisomerizations
	temp = (exp(-MeanPhotoisoms) * MeanPhotoisoms^n / factorial(n));
	% multiply by normalization factor for Gaussian
	temp = temp * ((2 * 3.141592 * (DarkNoiseSD^2 + n * beta(1)^2))^(-0.5));
	% multiply by Gaussian with mean n*single-amp and variance dark + n*single variance
	temp = temp .* exp(-((x - n * beta(2)).^2) / (2 * (DarkNoiseSD^2 + n * beta(1)^2)));
	% add to fit for all photon counts
	fit = fit + temp .* NumResponses .* xinc;
end
