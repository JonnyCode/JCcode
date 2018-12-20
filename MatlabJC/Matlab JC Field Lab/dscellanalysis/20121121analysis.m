% data001 BW-8-2-0.48-11111 7200 s
% 
% data002 Short DS 2000 s trigger-interval: 12 s   - 8*2*12*10
% (spatial-periods '(64))  
% (temporal-periods '(32 256))
% (directions '(0 45 90 135 180 225 270 315))




[datarun002] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/', 'data002-map/data002-map', 1, '/stimuli/s02', 0);
[datarun001] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/', 'data001-map/data001-map', 0, 0, 1);
cellids = intersect(datarun002.cell_ids, datarun001.cell_ids);

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);


%Put spatial receptive field, time course, EI, polar+raster plot, pulse
%plot of those cells in a figure and save each cell as one pdf file


path = '/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/data002-map/All Plots/';
addpath(path);
sp = 64;
temp = 256;
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);


for aa = 1:length(cellids)
cell_id = cellids(1, aa);    
    
h3 = figure;
%Receptive Field from data000
ax1 = subplot(3,5,[1 2], 'Parent', h3);
plot_rf(datarun001, cell_id, 'foa', ax1, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0])
set(ax1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)

%Time course from data000
ax4 = subplot(3,5,[6 7], 'Parent', h3);
plot_time_course(datarun001,cell_id, 'figure', ax4, 'clear', false); 
title ('Temporal Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax4,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax4, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%EI from data000
ax5 = subplot(3,5,[11 12], 'Parent', h3);
plot_ei(datarun001, cell_id, 'foa', ax5, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%Polar Plot & Spike Raster from data002
figure;
[T R Unew Vnew] = polar_plots_one(rho(:,spatPer),theta(:,spatPer),U(:,spatPer),V(:,spatPer),num,get_cell_indices(datarun002, cell_id));
close;
all = 1:length(datarun002.stimulus.combinations);
[an indic] = sortrows(StimComb, 3); %sort according to angle
all = all(indic);
[an indic] = sortrows(an, 2); %sort according to speed
all = all(indic);
[an indic] = sortrows(an, 1); %sort according to sp period
all = all(indic); %
subpin = [10 5 4 3 8 13 14 15]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,5,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 10, 'stop', 12, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,5,9)
%set(gca, 'FontName', 'Helvetica', 'FontSize', 11)
set(gca, 'FontSize', 9)
p = polar(T{tempPer,1}, R{tempPer,1});
h = findall(gca, 'type', 'line'); % find all lines in the current figure
h(h==p) = []; % remove the line you ploted from the list.
set(h, 'LineStyle', '--');
hold on;
h1 = compass(Unew(tempPer,1),Vnew(tempPer,1), 'r'); %Vector average plot
set(h1,'linewidth',2) 
hold off;
title ('Drifting Grating Rasters', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0], 'Position', [0.16, 5.3, 18.23])


save_figure_pdf(path, num2str(cellids(1,aa)), h3) %uncommnet

% pathname = ['/Users/sravi/matlab/DS cell analysis/2012-10-15-0/AllPlotsInOne-512/'  num2str(cellids(1,aa))];
% set(h3,'PaperOrientation','landscape');
% set(h3,'PaperUnits','normalized');
% set(h3,'PaperPosition', [0 0 1 1]);
% print(h3,'-dpdf', pathname);

close(h3); %uncomment
%set(gcf, 'Color', [1 1 1]);

end

chos = get_cell_indices(datarun002, cellids);

%Raster Plots for DG DS
all_rasterplots_pdf(NumSpikesCell,StimComb,datarun002, rho,theta,U,V,num, chos, 1,2);


plot(magin{1,1}, magin{2,1},'o','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);

ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cid = 1:length(X);
chos = cid(in);


%DS Cells
%chos = ds cells index wrt 324 cells
ds = cellids(1, chos); %cell ids of ds cells

NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, ds), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds] = dscellanalysis(NumSpikesCellDS, StimComb);

plot(angleds{1,1}, angleds{2,1}, '+')

ginput(1);
Y = get(get(gca,'Children'),'YData');
X = get(get(gca,'Children'),'XData');
h = impoly;
accepted_pos = wait(h);
hold on;
drawPolygon(accepted_pos(:,1), accepted_pos(:,2));
in = inpolygon(X, Y, accepted_pos(:,1), accepted_pos(:,2));
plot(X(in),Y(in),'r+');
cid = 1:length(X);
chos2 = cid(in);
ds(1,chos2)




A = [468 1115 1699 1896 2659 2747 3242 3437 4052 5463 6228]; %180 degrees ON
A = [2044  4384 4532 5252 5342 6619]; %180 degrees ON OFF
A = [468 1115 1699 1896 2659 2747 3242 3437 4052 5463 6228 2044  4384 4532 5252 5342 6619]; %180 degrees all
A = [1399 3603 6377 6979]; %-90 degrees
A = [1294 1956 2387 2642 3167 4037 4280 4730 5177 5387 5643]; %0 degrees
A = [107 380 736 904 920 1277 1368 1595 1788 2553 3425 3468 3812 4369 4411 4757 4922 5059 5225 5374 5881 5915 6274]; %90 degrees

plot_rf_summaries(datarun001, A,'label', true)
plot_rf_portraits(datarun001,A)
plot_time_courses(datarun001, A, 1, 1);
% datarun001 = get_autocorrelations(datarun001, 'all');
plot_autocorrelograms(datarun001, A, 'foa', 0);





















[tc nontc] = get_time_courses_matrix(datarun001, cellids);




[datarun002old] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/', 'data002-map-old/data002-map-old', 1, '/stimuli/s02', 0);
[datarun001old] = load_dsdata('/Analysis/sravi/Rat/Glaucoma/2012-11-21-1/data001-3600-7200s/','data001-map-old/data001-map-old', 0, 0, 1);

cellids = intersect(datarun002old.cell_ids, datarun001old.cell_ids);


[NumSpikesCell2, StimComb2] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
[NumSpikesCell2old, StimComb2old] = get_spikescellstim(datarun002old, datarun002old.cell_ids, 0);



[magin2 dsindexin2 magmaxin2 magavein2 anglein2 rhoin2 thetain2 numin2 Uin2 Vin2] = dscellanalysis(NumSpikesCell2, StimComb2);
[magin2old dsindexin2old magmaxin2old magavein2old anglein2old rhoin2old thetain2old numin2old Uin2old Vin2old] = dscellanalysis(NumSpikesCell2old, StimComb2old);


plot(magin2old{1,1}, magin2old{2,1},'o','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
