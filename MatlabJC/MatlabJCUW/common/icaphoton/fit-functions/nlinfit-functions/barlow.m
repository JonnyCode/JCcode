function [fit] = barlow(beta, x)

fit = abs(beta(1)) * ((abs(beta(2)) + x).^beta(3));
