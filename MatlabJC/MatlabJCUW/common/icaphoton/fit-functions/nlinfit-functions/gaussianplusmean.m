function fit = gaussianplusmean(beta, x)

fit = x - beta(2);
fit = beta(1) .* exp(-fit .* fit ./ (2 * (beta(3)^2))) / ((2 * 3.14159 * beta(3)^2)^(0.5)) + beta(4);
