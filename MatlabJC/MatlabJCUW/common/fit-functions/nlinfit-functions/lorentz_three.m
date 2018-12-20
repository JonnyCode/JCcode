function [fit] = lorentz_three(beta, x)

beta(2) = 5;
beta(4) = 30;
beta(6) = 535;

fit = abs(beta(1)) ./ (1 + (x ./ abs(beta(2))).^2).^2;
fit = fit + abs(beta(3)) ./ (1 + (x ./ abs(beta(4))).^2);
fit = (fit + abs(beta(5)) ./ (1 + (x ./ abs(beta(6))).^2));
