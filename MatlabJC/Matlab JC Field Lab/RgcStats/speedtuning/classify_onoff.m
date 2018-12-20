function [id_sub, idx_sub] = classify_onoff(mag_pca, ds_id, pc1, pc2, manual)

[id_sub, idx_sub] = deal(cell(2, 1));

FigHandle = figure;
set(FigHandle, 'Position', [1 1 380 400])

[~,scores,~,~] = princomp(mag_pca);
plot(scores(:, pc1), scores(:, pc2), 'o')
hold on
[x, y] = ginput;
plot(x, y)
IN = inpolygon(scores(:, pc1), scores(:, pc2), x, y);
[~, idx_init] = find(IN' == 1);
id_init = ds_id(idx_init);

id_sub{1} = setdiff(ds_id, id_init);
id_sub{2} = id_init;
idx_sub{1} = setdiff([1:length(ds_id)], idx_init);
idx_sub{2} = idx_init;

if ~manual
    [~, ~, ib] = intersect(id_init, ds_id);
    vc = ones(length(ds_id),1);
    vc(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

    X(:,1) = scores(:, pc1);
    X(:,2) = scores(:, pc2);
    [idx, ~, ~] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, vc);

    id_sub{1} = ds_id(idx==1);
    id_sub{2} = ds_id(idx==2);

    [~, idx_sub{1}] = intersect(ds_id, id_sub{1});
    [~, idx_sub{2}] = intersect(ds_id, id_sub{2});
end

end
