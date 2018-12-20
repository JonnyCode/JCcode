function [fit] = polynomial_2(beta, x)

fit = beta(1) + beta(2) .* x + beta(3) .* x^2;
