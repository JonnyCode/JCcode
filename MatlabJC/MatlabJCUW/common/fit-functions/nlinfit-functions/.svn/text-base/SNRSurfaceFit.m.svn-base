% fit SNR surface with separable model

function fit = SNRSurfaceFit(beta, x)

global MeanPhotoIsoms;
global MaxSNR;

yfit = normcdf(x, beta(1), abs(beta(2)));
xfit = normcdf(MeanPhotoIsoms, beta(3), abs(beta(4)));

for xpos = 1:length(xfit)
	for ypos = 1:length(yfit)
%		fit((xpos-1) * length(yfit) + ypos) = xfit(xpos) * yfit(ypos) * abs(beta(5));
		fit((xpos-1) * length(yfit) + ypos) = xfit(xpos) * yfit(ypos) * MaxSNR;
	end
end

fit = fit';
