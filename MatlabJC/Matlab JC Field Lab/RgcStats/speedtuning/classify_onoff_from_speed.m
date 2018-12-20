function [id_sub, idx_sub] = classify_onoff_from_speed(datarun, ds_id, pc1, pc2, manual)

% eg. manually classify with 1st and 2nd pc
% [id_sub, idx_sub] = classify_onoff(mag_pca, ds_id, 1, 2, 1)

% adapted by JC 2018-11-20-0 to point to XYs code versions 'Vxy'

[NumSpikesCell, MaxRate, StimComb]  = get_spikescellstimVxy(datarun,ds_id,0,1);
DG = sort_direction(dscellanalysisVxy(NumSpikesCell, StimComb, datarun));
MAG_all_norm_dg = normalize_MAG(DG);
[id_sub, idx_sub] = classify_onoff(MAG_all_norm_dg', ds_id, pc1, pc2, manual); % select the right cluster!!!

figure
plot(MAG_all_norm_dg(:,idx_sub{1}),'r')
hold on
plot(MAG_all_norm_dg(:,idx_sub{2}),'b')

end