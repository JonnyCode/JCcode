function [fit] = exponential(coef, x)

fit = coef(1) + coef(2) .* exp(-x ./ coef(3));
