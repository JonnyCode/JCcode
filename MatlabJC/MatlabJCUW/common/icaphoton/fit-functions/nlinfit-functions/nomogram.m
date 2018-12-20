function [fit] = nomogram(beta, x)

fit(1:length(x)) = 0;
temp = (561 * x / 1000) / beta(8);

for n=0:6
	fit = fit - beta(n+1) * ((log10(temp)').^n);
end

fit = fit';
