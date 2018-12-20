function [fit] = poisson_filter_sine(beta, x)

fit = beta(1) .* (x/abs(beta(3))).^abs(beta(2)) .* exp(-x / abs(beta(3))) .* sin(2 * 3.14159 * x / abs(beta(4)));
