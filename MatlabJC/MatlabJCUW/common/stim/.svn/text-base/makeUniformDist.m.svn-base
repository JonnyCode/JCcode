function M = makeUniformDist(M,C)
%C is contrast
M_flat = M(:);
bins = [0 prctile(M_flat,[1:1:100])];
M_orig = M;
for i=1:length(bins)-1
    M(M_orig>bins(i) & M_orig<=bins(i+1)) = i*(C/100);
end
M = M + .5 - C/2;
%M = M ./ (mean(M(:))./.5);