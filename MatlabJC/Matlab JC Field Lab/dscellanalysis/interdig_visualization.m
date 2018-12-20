% compute and plot contours for a mosaic of RFs

data_spec = '/Volumes/lab/Experiments/Array/Analysis/Chichilnisky-lab/2012-05-31-1/sr-map/data001-map/data001-map';
cell_spec = {22};



subplot(1,2,1)
imagesc(datarun000.stas.rfs{82,1})
subplot(1,2,2)
imagesc(datarun000.stas.summaries_filt{1,82})

% initialize struct
datarun000 = load_data(data_spec);

% load STAs
datarun000 = load_sta(datarun000,'load_sta',[],'keep_java_sta',true);

% load cell types
datarun000 = load_params(datarun000,'verbose',1);

% compute RFs, etc
datarun000 = get_sta_summaries(datarun000,cell_spec);

% filter summary frames (saved to datarun000.stas.summaries_filt{})
datarun000 = get_rfs_filtered(datarun000,cell_spec,'verbose',1,'filt_params',struct('filt_type','gauss','radius',.4));

% Do Delaunay Triangulation; for determining RoI for calculating UI
datarun000 = do_delaunay_tri(datarun000, cell_spec, true);
    

% Cut triangles with excessively long edges; indicates where cells are
% likely missing in mosaic.
datarun000 = cull_delaunay_tri(datarun000, cell_spec, 1.9, 'plot', true);

% Build Region of Interest from culled Delaunay Triangulation
datarun000 = build_rf_roi(datarun000, cell_spec,1);

% Get contour threshold that maximizes uniformity index
datarun000 = maximize_uniformity_index(datarun000, cell_spec, 0.3, 0.05, 0.01);

% compute contours (summary frames normalized, then saved to datarun000.stas.contours)
datarun000 = get_rf_contours(datarun000, cell_spec, datarun000.stas.mosaics{cell_spec{:}}.best_thresh, 'norm_params', struct('method', 'peak'), 'verbose', 1);


% plot contours
plot_rf_summaries(datarun000, cell_spec, 'plot_contours', 1, 'foa', 1, 'label', 0);
lock;

% plot simplified contours, with alpha fills
datarun000 = simplify_rf_contours(datarun000, cell_spec);
plot_rf_summaries(datarun000, cell_spec, 'plot_contours', 1, 'contours_field', 'rf_contours_simple', 'foa', 1, 'label', 0, 'contour_fill', 'r', 'contour_alpha', 0.5);



% ALTERNATIVELY...
figure
plot_rf_coloring(datarun000, cell_spec, 'rfs', 'summaries_filt');
plot_rf_coloring(datarun000, cell_spec);

%2012-10-31

ont2cellind = get_cell_indices(datarun000, ont2);
for i = 1:length(ont2)
    imagesc(datarun000.stas.rfs{ont2cellind(i),1});
    colormap(gray);
    hold on;
    plot(datarun000.stas.rf_contours_simple{1,ont2cellind(i)}{1,1}.x, datarun000.stas.rf_contours_simple{1,ont2cellind(i)}{1,1}.y, 'r','LineWidth',2);
    pause;
    close;
end

ont3cellind = get_cell_indices(datarun000, ont3);
for i = 1:length(ont3)
    imagesc(datarun000.stas.rfs{ont3cellind(i),1});
    colormap(gray);
    hold on;
    plot(datarun000.stas.rf_contours_simple{1,ont3cellind(i)}{1,1}.x, datarun000.stas.rf_contours_simple{1,ont3cellind(i)}{1,1}.y, 'r','LineWidth',2);
    pause;
    close;
end

ont2cellind = get_cell_indices(datarun000, [ont1] );
for i = 1:length(ont2cellind)
    hold on;
    plot(datarun000.stas.rf_contours_simple{1,ont2cellind(i)}{1,1}.x, datarun000.stas.rf_contours_simple{1,ont2cellind(i)}{1,1}.y, 'r','LineWidth',1.2);
    pause;
end

ont3cellind = get_cell_indices(datarun000, ont3);
for i = 1:length(ont3)
    hold on;
    plot(datarun000.stas.rf_contours_simple{1,ont3cellind(i)}{1,1}.x, datarun000.stas.rf_contours_simple{1,ont3cellind(i)}{1,1}.y, 'r','LineWidth',1.2);
    title(num2str(ont3cellind(i)))
    pause;
end
    
ont2extra = [1606 4113 4156 4786 6888];
ont2extra = [ont2 2086 2371 7669 4186];

ont2extra = [454 496 528 1578 2716 3152 3530 4145 4501 6155];
ont2extracellind = get_cell_indices(datarun000, nott1);
for i = 1:length(ont2extracellind)
    hold on;
    plot(datarun000.stas.rf_contours_simple{1,ont2extracellind(i)}{1,1}.x, datarun000.stas.rf_contours_simple{1,ont2extracellind(i)}{1,1}.y, 'b','LineWidth',1.2);
    pause;
end

datarun000 = cull_delaunay_tri(datarun000, cell_spec, 2.7, 'plot', true);
hold on;
plot_rf_summaries(datarun000, cell_spec, 'plot_contours', 1, 'contours_field', 'rf_contours_simple', 'foa', 1, 'label', 0, 'contour_fill', 'r', 'contour_alpha', 0.5);

%2012-10-15
plot_rf_summaries(datarun000, cell_spec,'label', true)
plot_time_courses(datarun000, [ont1 [1969	3721	4096	4427	5104	5941 4501 3152 3530]], 'all' , true, 'bw', true, 'normalize', true);
   

