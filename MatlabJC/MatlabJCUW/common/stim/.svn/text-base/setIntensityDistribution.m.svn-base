function M = setIntensityDistribution(M,D)
M_flat = M(:);
nBins = length(D);
minBin = 100/nBins;
bins = [0 prctile(M_flat,linspace(minBin,100,nBins))];
M_orig = M;
for i=1:length(bins)-1
    M(M_orig>bins(i) & M_orig<=bins(i+1)) = D(i);
end
%M = M ./ (mean(M(:))./.5);