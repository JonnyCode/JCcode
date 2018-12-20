function D = getIntensityDistribution(M,nBins)
M_flat = M(:);
minBin = 100/nBins;
bins = [0 prctile(M_flat,linspace(minBin,100,nBins))];
for i=1:length(bins)-1
    D(i) = mean(mean(M(M>bins(i) & M<=bins(i+1))));
end

D = naninterp(D);