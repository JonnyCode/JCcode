function [fit] = exponential(coef, x)

fit = 1 - coef(1) .* log(x .* coef(2)) ./ (1 + x ./ (coef(3) * 1e4));
