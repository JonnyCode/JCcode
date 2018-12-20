% fit SNR surface with separable model

function fit = SNRSurfaceFit2(beta, x)

global MeanPhotoIsoms;
global MaxSNR;

yfit = normcdf(x, beta(1), abs(beta(2)));
xfit = normcdf(MeanPhotoIsoms, beta(3), abs(beta(4)));

fit = MaxSNR * xfit .* yfit;
%fit = beta(5) * xfit .* yfit;

fit = fit';
