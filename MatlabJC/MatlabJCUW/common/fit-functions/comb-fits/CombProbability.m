function [ProbabilityA, ProbabilityB] = CombProbability(CombParameters, Resp, RespInc)
% CombProbability
%
% Probability of response Resp for comb distributions A and B from CombParameters.

MaxHit = 10;

ProbA = 0;
ProbB = 0;
temp = Resp;

for n = 0:MaxHit
	% Poisson factor for probability of n photoisomerizations
	temp = (exp(-CombParameters.RateA) * CombParameters.RateA^n / factorial(n));
	% multiply by normalization factor for Gaussian
	temp = temp * ((2 * 3.141592 * (CombParameters.DarkNoiseSD^2 + n * CombParameters.UnitarySDA^2))^(-0.5));
	% multiply by Gaussian with mean n*single-amp and variance dark + n*single variance
	temp = temp .* exp(-((Resp - n * CombParameters.UnitaryAmpA).^2) / (2 * (CombParameters.DarkNoiseSD^2 + n * CombParameters.UnitarySDA^2)));
	% add to fit for all photon counts
	ProbA = ProbA + temp * RespInc;
	% Poisson factor for probability of n photoisomerizations
	temp = (exp(-CombParameters.RateB) * CombParameters.RateB^n / factorial(n));
	% multiply by normalization factor for Gaussian
	temp = temp * ((2 * 3.141592 * (CombParameters.DarkNoiseSD^2 + n * CombParameters.UnitarySDB^2))^(-0.5));
	% multiply by Gaussian with mean n*single-amp and variance dark + n*single variance
	temp = temp .* exp(-((Resp - n * CombParameters.UnitaryAmpB).^2) / (2 * (CombParameters.DarkNoiseSD^2 + n * CombParameters.UnitarySDB^2)));
	% add to fit for all photon counts
	ProbB = ProbB + temp * RespInc;
end

ProbabilityA = ProbA;
ProbabilityB = ProbB;
