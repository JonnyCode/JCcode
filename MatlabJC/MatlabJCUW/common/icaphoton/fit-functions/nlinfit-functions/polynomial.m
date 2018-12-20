function [fit] = polynomial(beta, x)

NumCoef = length(beta);

for coef = 1:NumCoef
	if (coef == 1)
		fit = beta(coef) * x;
	else
		fit = fit + beta(coef) * x.^coef;
	end
end
fit = fit + 1;

	
