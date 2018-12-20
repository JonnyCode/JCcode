function[dsindex] = DSI_finder(mag, rho)

% Function calculates and returns direction selective index, dsi =
% Vector_sum/ sum of vector mags

% modified from SR code "DS_index_one"
% JC
% Last revision: 3/13/15

dsindex = cell(length(mag),1);
for i = 1:length(mag)
    dsindex{i,1} = mag{i,1}./sum(rho{i,1},2)';
    [dsindex{i,1} ntimp] = exciseRows(dsindex{i,1}', zeros(size(dsindex{i,1},2),1));
end
end