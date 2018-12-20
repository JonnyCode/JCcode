function [fit] = poisson_filter(beta, x)

x1 = x(:,1);
fit = beta(1) .* (x1/beta(2)).^beta(3) .* exp(-x1 / beta(2));
