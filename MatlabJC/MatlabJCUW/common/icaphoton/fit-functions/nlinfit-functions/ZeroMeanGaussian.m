function fit = ZeroMeanGaussian(beta, x)

global NumResponses;
xinc = x(2) - x(1);

fit = x;
fit = NumResponses .* xinc .* exp(-fit .* fit ./ (2 * (beta(1)^2))) / (2 * 3.14159 * beta(1)^2)^(0.5);
