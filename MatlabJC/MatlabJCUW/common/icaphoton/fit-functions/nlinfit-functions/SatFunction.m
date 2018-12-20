function fit = SatFunction(beta,x)

%	fit = SatFunction(beta,x)
%
%	saturating function for determining the shape of the surface

fit = beta(1) .* (1 - exp(-(beta(2) .* x))) + 0.5;
