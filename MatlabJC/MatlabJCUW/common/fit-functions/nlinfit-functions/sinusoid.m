function [fit] = sinusoid(coef, x)

fit = coef(1) .* sin(2 * 3.14159 * x * 0.001 + coef(3));
