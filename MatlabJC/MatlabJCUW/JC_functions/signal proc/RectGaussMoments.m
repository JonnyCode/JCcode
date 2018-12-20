function [RectGaussMean, RectGaussSD] = RectGaussMoments(GaussMean, GaussSD, RectThresh, SampStep)

Samples = [RectThresh:SampStep:GaussMean+5*GaussSD];

NormFact = sum(normpdf(Samples, GaussMean, GaussSD));

RectGaussMean = sum(Samples .* normpdf(Samples, GaussMean, GaussSD)) / NormFact;
RectGaussVar = sum(Samples.^2 .* normpdf(Samples, GaussMean, GaussSD)) / NormFact;
RectGaussSD = sqrt(RectGaussVar - RectGaussMean^2);
