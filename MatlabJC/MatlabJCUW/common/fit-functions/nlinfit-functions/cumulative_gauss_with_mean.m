function [fit] = cumulative_gauss_with_mean(beta, x)

fit = normcdf(x, beta(2), abs(beta(1)));
