function fit = normgaussian(beta, x)

fit = x - beta(1);
fit = exp(-fit .* fit ./ (2 * (beta(2)^2))) / ((2 * 3.14159 * beta(2)^2)^(0.5));

if (sum(fit) > 0)
    fit = fit / sum(fit);
end
fit = fit * 0.5 * 1000;
