function [fit] = lorentzsum(beta, x)

fit = abs(beta(1)) ./ (1 + (x ./ abs(beta(2))).^2).^2;
fit = (fit + abs(beta(3)) ./ (1 + (x ./ abs(beta(4))).^2));
