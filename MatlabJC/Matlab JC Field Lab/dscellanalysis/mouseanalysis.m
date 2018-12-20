[datarun000] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data000-3600-7200-dn-sr-map/data000-3600-7200-dn-sr-map', 0, 0, 1);
[datarun002] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data004-map/data004-map', 1, '/stimuli/s04', 0);
[datarun003] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data002-map/data002-map', 0, 0, 0);
[datarun001] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data001-60-600-dn-sr-map/data001-60-600-dn-sr-map', 0, 0, 0);

cellids = intersect(intersect(intersect(datarun000.cell_ids, datarun001.cell_ids), datarun002.cell_ids), datarun003.cell_ids);
cellids = intersect(datarun000.cell_ids, datarun002.cell_ids);

c1 = intersect(intersect(datarun000.cell_ids, datarun001.cell_ids), datarun002.cell_ids);
c2 = intersect(datarun000.cell_ids, datarun002.cell_ids);

[datarun0000409] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data021-map/data021-map', 0, 0, 1);


allplotsinone(datarun000, datarun001, datarun002, cellids1, rho,theta,U,V,num,StimComb, 64, 32, '/Analysis/sravi/Mouse/2013-06-05-0/AllPlots/', wh, gr);

[datarun003] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data003-map/data003-map', 1, '/stimuli/s03', 0);
[datarun001] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data001-map/data001-map', 0, 0, 0);

[datarun004560] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data004-560-3360-dn-sr-map/data004-560-3360-dn-sr-map', 1, '/stimuli/s04', 0);



[datarun300] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data000-300-7200-dn-sr-map/data000-300-7200-dn-sr-map', 0, 0, 1);

[datarun600] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data000-600-7200-dn-sr-map/data000-600-7200-dn-sr-map', 0, 0, 1);

[datarun1800] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data000-1800-7200-dn-sr-map/data000-1800-7200-dn-sr-map', 0, 0, 1);

ntersect(intersect(datarun000.cell_ids, datarun300.cell_ids), datarun001.cell_ids);

cellids = intersect(datarun000.cell_ids, datarun004.cell_ids);
intersect(datarun1800.cell_ids, datarun004.cell_ids)

intersect(datarun3600.cell_ids, datarun001.cell_ids)
intersect(datarun1800.cell_ids, datarun001.cell_ids)

intersect(datarun3600.cell_ids, datarun002.cell_ids)
intersect(datarun1800.cell_ids, datarun002.cell_ids)



intersect(datarun3600.cell_ids, datarun003.cell_ids)
intersect(datarun1800.cell_ids, datarun003.cell_ids)



intersect(datarun3600.cell_ids, datarun600.cell_ids)
intersect(datarun3600.cell_ids, datarun1800.cell_ids)

intersect(datarun3600.cell_ids, datarun004.cell_ids)
intersect(datarun3600.cell_ids, datarun004560.cell_ids)


[datarun004560] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data004-560-3360-dn-sr-map/data004-560-3360-dn-sr-map', 1, '/stimuli/s04', 0);


[datarun00260] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data002-60-600-dn-sr-map/data002-60-600-dn-sr-map', 0, 0, 0);

intersect(datarun3600.cell_ids, datarun002.cell_ids)
intersect(datarun3600.cell_ids, datarun00260.cell_ids)

intersect(datarun3600.cell_ids, datarun004560.cell_ids)

cellids = intersect(intersect(intersect(datarun3600.cell_ids, datarun00160.cell_ids), datarun004.cell_ids), datarun002.cell_ids);
cellids = intersect(intersect(datarun000.cell_ids, datarun300.cell_ids), datarun001.cell_ids);
cellids = intersect(intersect(intersect(datarun000.cell_ids, datarun001.cell_ids), datarun002.cell_ids), datarun003.cell_ids);

[NumSpikesCell, StimComb] = get_spikescellstim(datarun002, datarun002.cell_ids, 0);
NumSpikesCellIn = NumSpikesCell(get_cell_indices(datarun002, cellids), :);
[mag dsindex magmax magave angle rho theta num U V] = dscellanalysis(NumSpikesCell, StimComb);
[magin dsindexin magmaxin magavein anglein rhoin thetain numin Uin Vin] = dscellanalysis(NumSpikesCellIn, StimComb);

all_rasterplots_pdf(NumSpikesCell,StimComb,datarun004, rho,theta,U,V,num, get_cell_indices(datarun004, cellids), 1,2)

a = 1;

for i = 1:length(numin)
    for j = i+1:length(numin)
 subplot(3,7, a)
 plot(magin{i,1}, magin{j,1}, 'o')
 xlabel(num2str(numin(i,1)));
 ylabel(num2str(numin(j,1)));
 a = a+1;
    end
end



a = 1;

for i = 1:length(numin)
    for j = i+1:length(numin)
 plot(magin{i,1}, magin{j,1}, 'o')
 xlabel(num2str(numin(i,1)));
 ylabel(num2str(numin(j,1)));
 a = a+1;
 pause;
    end
end

plot(magin{2,1}, magin{3,1}, 'o')

plot(magin{2,1}, magin{3,1},'o','Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 7, 'MarkerEdgeColor' , [0 0 1], 'MarkerFaceColor' , [.7 0.8 1]);
title ('Normalized vector average in preferred direction', 'FontName','Tahoma', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
xlabel('Fast temporal period - 16', 'FontName','AvantGarde', 'FontSize', 16, 'Color', [0 0 0]);
ylabel('Slow temporal period - 32', 'FontName','AvantGarde', 'FontSize', 16, 'Color', [0 0 0]);
set(gca, 'FontName', 'Helvetica', 'FontSize', 16);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [.005 .005],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'DS', gcf)



plot(magin{2,1}, magin{3,1}, 'o')
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


ds = cellids(1, chos); %cell ids of ds cells

NumSpikesCellDS = NumSpikesCell(get_cell_indices(datarun002, ds), :);
[magds dsindexds magmaxds magaveds angleds rhods thetads numds Uds Vds] = dscellanalysis(NumSpikesCellDS, StimComb);

polar_plots_all(Uds,Vds,numds, magds)
hh = title ('Preferred Directions of cells', 'FontName','Tahoma', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(hh,'Position',[6,8,1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'Polar', gcf)

plot(angleds{2,1}, angleds{3,1}, '+');
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
deg = ds(1,chos);

%90_degs = [317 662 1262 3976 4490 4805 4847 5271];
%0_degs = [422 514 2417 2433 2451 2659 2944 3184 3214 3602 3754 4802 5361 5671 5749 6106 7533 7593];
%180_degs = [181 617 647 1517 1636 1790 2297 2342 2358 2462 2688 2720 3142 3601 3919 4114 4278 4804 4834 5088 5508 6125 6316 6484 7322];
%270_degs = [963 1143 1639 1803 2180 2613 3377 3422 3766 6421 6513 6963 7143 7293];

A = [963 1143 1639 1803 2180 2613 3377 3422 3766 6421 6513 6963 7143 7293];
plot_rf_summaries(datarun000, A, 'coordinates', 'monitor');
title ('Receptive Field Mosaic', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0])
xlabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);
ylabel('pixels', 'FontName', 'AvantGarde', 'FontSize', 18,'Color', [0 0 0]);  
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gcf, 'Color', [1 1 1]);
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'mosaic270deg', gcf)





xlabel('16')
ylabel('32')
save_figure_pdf('/Analysis/sravi/Mouse/2013-06-05-0/DS Plots/', 'data004-best-classf',gcf)

save_figure_pdf('/Analysis/sravi/Mouse/2013-06-05-0/DS Plots/', 'data004-classf',gcf)











subplot(1,2,1)
plot(magin{1,1}, magin{2,1}, 'o');
xlabel('TP 32');
ylabel('TP 64');
title('SP 64');
subplot(1,2,2)
plot(magin{1,2}, magin{2,2}, 'o');
xlabel('TP 32');
ylabel('TP 64');
title('SP 128');
save_figure_pdf('/Analysis/sravi/Mouse/2013-06-05-0/DS Plots/', 'data003-classf',gcf)


plot(magin{2,1}, magin{3,1}, 'o');

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



sp = 64;
temp = 32;
path = '/Analysis/sravi/Mouse/2013-06-05-0/AllPlots/';
spatPer = find(unique(StimComb(:,1))==sp);
tempPer = find(unique(StimComb(:,2))==temp);



for aa = 1:length(cellids)
cell_id = cellids(1, aa);    
    
h3 = figure;
%Receptive Field from data000
ax1 = subplot(3,9,[1 2], 'Parent', h3);
plot_rf(datarun000, cell_id, 'foa', ax1, 'scale', 5);
title ('Spatial Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0])
set(ax1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)

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
subpin = [14 5 4 3 12 21 22 23]; %Places on subplot for each angle from 0 to 315 deg
in1 = ismember(an(:,1), sp)'; % spatial period
SC2 = an(ismember(an(:,1), sp), :);% spatial period
A1 = all(in1);%

in = ismember(SC2(:,2), num(tempPer))'; %temporal period
A = A1(in);
for  j = 1:sum(in)
    trigpre = ismember(datarun002.stimulus.trial_list,A(j));
    destaxes=subplot(3,9,subpin(1,j),'Parent', h3);
    spikesbytrials = get_raster(datarun002.spikes{get_cell_indices(datarun002, cell_id),1}, datarun002.stimulus.triggers(trigpre), 0, 8, 0, 6, 'stop', 10, 'foa', destaxes, 'tic_color', [0 0 0]);
end
subplot(3,9,13)
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


%Pulse Raster from data001
ax3 = subplot(3,9,[6 7 15 16 24 25], 'Parent', h3);
wh = datarun001.triggers(21:4:length(datarun001.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr = datarun001.triggers(22:4:length(datarun001.triggers),1); %Might change with dataset
[hp ap] = pulse_analysis(datarun001, get_cell_indices(datarun001, cell_id), 0, path, wh, gr); %change raster plot axes
fig3 = get(ap,'children');
copyobj(fig3,ax3);
close(hp);
title ('W-G-B-G', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]); 
%set(ax3, 'FontName', 'Helvetica', 'FontSize', 12);
set(ax3,'Box', 'off', 'TickDir','out', 'TickLength', [.005 .005],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
%ylabel('trials', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%Pulse Raster from data003
ax3 = subplot(3,9,[8 9 17 18 26 27], 'Parent', h3);
wh3 = datarun003.triggers(1:4:length(datarun003.triggers), 1); %For 2012-10-10-1 dataset, that is how the triggers are arranges
gr3 = datarun003.triggers(2:4:length(datarun003.triggers),1); %Might change with dataset
[hp ap] = pulse_analysis(datarun003, get_cell_indices(datarun003, cell_id), 0, path, wh3, gr3); %change raster plot axes
fig3 = get(ap,'children');
copyobj(fig3,ax3);
close(hp);
title ('W-G-DG-G', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]); 
%set(ax3, 'FontName', 'Helvetica', 'FontSize', 12);
set(ax3,'Box', 'off', 'TickDir','out', 'TickLength', [.005 .005],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
%ylabel('trials', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%Time course from data000
ax4 = subplot(3,9,[10 11], 'Parent', h3);
plot_time_course(datarun000,cell_id, 'figure', ax4, 'clear', false); %'clear', false needed at varargin param?
title ('Temporal Receptive Field', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax4,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax4, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


%ACF from data000
ax5 = subplot(3,9,[19 20], 'Parent', h3);
plot(datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.bins, datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.probabilities);
title ('ACF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


save_figure_pdf(path, num2str(cellids(1,aa)), h3) %uncommnet

close(h3); %uncomment

end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
%Plotting cells with only WN and DG


for aa = 1:length(cellids)
cell_id = cellids(1, aa);    
    
h3 = figure;
%Receptive Field from data000
ax1 = subplot(2,2,1, 'Parent', h3);
plot_rf(datarun000, cell_id, 'foa', ax1, 'scale', 5);
title ('Spatial RF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0])
set(ax1,'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)


%Time course from data000
ax4 = subplot(2,2,3, 'Parent', h3);
plot_time_course(datarun000,cell_id, 'figure', ax4, 'clear', false); %'clear', false needed at varargin param?
title ('Temporal RF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax4,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax4, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('frame number', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);

%ACF from data000
ax5 = subplot(2,2,2, 'Parent', h3);
plot(datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.bins, datarun000.autocorrelation{get_cell_indices(datarun000, cell_id),1}.probabilities);
title ('ACF', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax5,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax5, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('bins', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('probability', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);



%EI from data000
ax6 = subplot(2,2,4, 'Parent', h3);
plot_ei(datarun000, cell_id, 'foa', ax6, 'coordinates' , 'array');
title ('Electrophysiological Image', 'FontName','Tahoma', 'FontSize',12,'FontWeight','bold', 'Color', [0 0 0]);
set(ax6,'Box', 'off', 'TickDir','out', 'TickLength', [.01 .01],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(ax6, 'FontName', 'Helvetica', 'FontSize', 12);
xlabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);
ylabel('microns', 'FontName','AvantGarde', 'FontSize', 11, 'Color', [0 0 0]);


save_figure_pdf(path, num2str(cellids(1,aa)), h3) %uncommnet

close(h3); %uncomment

end


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 2013-06-05-0 cells of each of the 5off and 2 on types
 dir('/Analysis/sravi/Mouse/2013-06-05-0/AllPlots/OFF/t1')

ca = [2000 4306 6469 7532 2643 3483 4698 6872 3031 3707 5191 6888 1507 3201 4171 5266 6935 1730 3287 6138 7459];  %on t1     
cb = [1982 286 4492 6243 7277 2086 3093 4741 6391 7682 1249 2161 3152 512 902 1309 3543 5192 6782 948 1653 2611 4277 5927 6962 1834 2672 4339 6107 7069]; %on t2



ccold = [1727 2405 3676 4967 7292 1757 2506 3736 5043 7637 1083 1833 2644 3767 5268 1141 2030 2822 4246 5404 1172 2102 2958 4322 5566 1306 2177 3469 4336 5899 1699 2346 3587 4893 6605];      


cc = [1172 1699 1727 2405 3587 3676 3736 5404 5899 7292]; %off t1 further divided from pulses
cd = [1457 3168 436  5509 7036 1966 3376 4531 604 7264 1068 2071 3408 4879 6047 766 1097 2253 3691 5284 6198 7711 1171 303 6395]; %off t2
ce = [2716 3827 5221 1758 31 3994 5269 6632 1382 2132 3122 4022 527 6857 2284 3124 4097 5402 7186 1503 2537 3436 4471 557 7235 1651 2582 3513 4711 603 7291 1742 2629 3680 5146 6376];%off t3   
cf = [376 3631 4141 4426 4517 4756 4846 4906 5161 5327 5552 6587 7157 7504 7654]; %offt4
cg = [17 484 783 1351 1952 3063 3812 3946 4051 4669 6016 6796 7729]; %offt5

off = [17 31 108 122 303 316 376 436 453 466 484 527 557 603 604 766 783 827 872 886 904 905 931 980 991 1006 1068 1083 1084 1097 1141 1144 1171 1172 1201 1202 1246 1306 1351 1367 1382 1413 1426 1444 1456 1457 1471 1503 1546 1550 1564 1579 1651 1699 1727 1742 1756 1757 1758 1771 1833 1861 1877 1952 1966 2030 2071 2102 2118 2132 2177 2238 2253 2284 2326 2343 2346 2386 2405 2506 2521 2537 2582 2629 2644 2674 2701 2716 2717 2822 2882 2911 2958 3002 3063 3076 3122 3124 3151 3168 3211 3227 3242 3257 3276 3305 3331 3349 3376 3408 3436 3469 3470 3513 3515 3529 3587 3631 3633 3661 3676 3680 3691 3736 3767 3812 3826 3827 3946 3994 4021 4022 4051 4053 4071 4097 4141 4174 4203 4246 4248 4308 4322 4325 4336 4426 4428 4471 4486 4516 4517 4531 4669 4711 4712 4756 4817 4846 4864 4879 4893 4906 4967 5043 5087 5104 5106 5146 5161 5177 5221 5268 5269 5284 5327 5341 5356 5358 5402 5404 5432 5448 5477 5509 5552 5566 5570 5656 5672 5688 5702 5899 5911 6016 6047 6126 6152 6167 6198 6301 6332 6376 6395 6559 6587 6604 6605 6632 6721 6796 6857 6917 7036 7157 7186 7204 7235 7247 7264 7291 7292 7351 7504 7563 7637 7651 7654 7697 7711 7729];
ds = [1639 2342 2688 3377 3919 4802 5271 5750 6484 7096 1786 2358 2720 3422 3976 4804 5360 6079 6513 7143 1143 1790 2401 2941 351 4068 4805 5361 6106 662 7293 1262 1803 2417 2944 3544 4114 4834 5447 6125 6829 7322 1278 181 2433 3050 3601 422 4847 5462 617 6873 7533 1352 1816 2451 3142 3602 4278 4968 5508 6316 6919 7593 1517 1983 2462 317 3754 4487 5088 5643 6347 692 7657 1533 2180 2613 3184 3766 4488 5089 5671 6421 6963 77 1636 2297 2659 3214 3842 4490 514 5749 647 6976 963];  
on = [ 2000 4306 6469 7532 2643 3483 4698 6872 3031 3707 5191 6888 1507 3201 4171 5266 6935 1730 3287 6138 7459 1982  286 4492 6243 7277 2086 3093 4741 6391 7682 1249 2161 3152  512  902 1309 3543 5192 6782  948 1653 2611 4277 5927 6962 1834 2672 4339 6107 7069];
offclassf = [ 1172 1699 1727 2405 3587 3676 3736 5404 5899 7292 1457 3168 436 5509 7036 1966 3376 4531  604 7264 1068 2071 3408 4879 6047  766 1097 2253 3691 5284 6198 7711 1171  303 6395 2716 3827 5221 1758   31 3994 5269 6632 1382 2132 3122 4022  527 6857 2284 3124 4097 5402 7186 1503 2537 3436 4471  557 7235 1651 2582 3513 4711  603 7291 1742 2629 3680 5146 6376  376 3631 4141 4426 4517 4756 4846 4906 5161 5327 5552 6587 7157 7504 7654 17 484  783 1351 1952 3063 3812 3946 4051 4669 6016 6796 7729];







%%%%%%%%%%%%%%%%%%%%%%%%
 
 
[datarun0000605] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data000-3600-7200-dn-sr-map/data000-3600-7200-dn-sr-map', 0, 0, 1);
[datarun0000331] = load_dsdata('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/', 'data001-map/data001-map', 0, 0, 1);
[datarun0000409] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data021-map/data021-map', 0, 0, 1);


[datarun0010605] = load_dsdata('/Analysis/sravi/Mouse/2013-06-05-0/', 'data001-60-600-dn-sr-map/data001-60-600-dn-sr-map', 0, 0, 0);
[datarun0010331] = load_dsdata('/Analysis/sravi/Mouse/2013-03-31-0/data001/data001-2700-5341s/', 'data002-map/data002-map', 0, 0, 0);
[datarun0010409] = load_dsdata('/Analysis/sravi/Mouse/2013-04-09-0/', 'data020-map/data020-map', 0, 0, 0);



[tc0331 nontc] = get_time_courses_matrix(datarun0000331, datarun0000331.cell_ids);
[tc0409 nontc] = get_time_courses_matrix(datarun0000409, datarun0000409.cell_ids);
[tc0605 nontc] = get_time_courses_matrix(datarun0000605, datarun0000605.cell_ids);
x = 1:size(datarun0000409.stas.time_courses{1}, 1);

  time_course = time_course./norm(time_course); %normalize time course
  
  plot(x,tc0331(:,get_cell_indices(datarun0000331, 4997))/norm(tc0331(:,get_cell_indices(datarun0000331, 4997))));
  hold all;
  plot(x,tc0409(:,get_cell_indices(datarun0000409, 2207))/norm(tc0409(:,get_cell_indices(datarun0000409, 2207))));
  plot(x,tc0605(:,get_cell_indices(datarun0000605, 3587))/norm(tc0605(:,get_cell_indices(datarun0000605, 3587))));
   h = findall(gca, 'type', 'line');
  set(h, 'LineWidth', 2);
title ('Temporal Receptive Fields of OFF T2', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time to spike (ms)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('contrast', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
  legend('0331', '0409', '0605')
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'offt2-same-tc', gcf)

datarun0000331 = get_autocorrelations(datarun0000331, 'all');
datarun0000409 = get_autocorrelations(datarun0000409, 'all');
datarun0000605 = get_autocorrelations(datarun0000605, 'all');

plot(datarun0000331.autocorrelation{get_cell_indices(datarun0000331, 4997),1}.bins, datarun0000331.autocorrelation{get_cell_indices(datarun0000331, 4997),1}.probabilities);
hold all;
plot(datarun0000409.autocorrelation{get_cell_indices(datarun0000409, 2207),1}.bins, datarun0000409.autocorrelation{get_cell_indices(datarun0000409, 2207),1}.probabilities);
plot(datarun0000605.autocorrelation{get_cell_indices(datarun0000605, 3587),1}.bins, datarun0000605.autocorrelation{get_cell_indices(datarun0000605, 3587),1}.probabilities);

h = findall(gca, 'type', 'line');
set(h, 'LineWidth', 2);
title ('Autocorrelation Functions of OFF T2', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]);
set(gca,'Box', 'off', 'TickDir','out', 'TickLength', [0.02 0.02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
set(gca, 'FontName', 'AvantGarde', 'FontSize', 18);
xlabel('Time Difference (s)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('# of pairs', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
legend('0331', '0409', '0605')
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'offt2-same-acf', gcf)




[h0331, a0331, spikesbytrials0331, sumSpTrTrig0331, nhist0331] = pulse_analysis(datarun0010331, get_cell_indices(datarun0010331, 4997), 0, path, datarun0010331.triggers(1:4:length(datarun0010331.triggers), 1), datarun0010331.triggers(2:4:length(datarun0010331.triggers), 1)); %change raster plot axes
psth0331 = sum(nhist0331{1,1})/size(nhist0331{1,1},1);
psth0331 = psth0331./0.1;

[h0409, a0409, spikesbytrials0409, sumSpTrTrig0409, nhist0409] = pulse_analysis(datarun0010409, get_cell_indices(datarun0010409, 2207), 0, path, datarun0010409.triggers(1:4:length(datarun0010409.triggers), 1), datarun0010409.triggers(2:4:length(datarun0010409.triggers), 1)); %change raster plot axes
psth0409 = sum(nhist0409{1,1})/size(nhist0409{1,1},1);
psth0409 = psth0409./0.1;


[h0605, a0605, spikesbytrials0605, sumSpTrTrig0605, nhist0605] = pulse_analysis(datarun0010605, get_cell_indices(datarun0010605, 3587), 0, path, datarun0010605.triggers(21:4:length(datarun0010605.triggers), 1), datarun0010605.triggers(22:4:length(datarun0010605.triggers), 1)); %change raster plot axes
psth0605 = sum(nhist0605{1,1})/size(nhist0605{1,1},1);
psth0605 = psth0605./0.1;


stairs(0.1:0.1:12,psth0331);
hold all;
stairs(0.1:0.1:12,psth0409);
stairs(0.1:0.1:12,psth0605);

 h = findall(gca, 'type', 'line');
 set(h, 'LineWidth', 2);
title ('Full Field Pulse PSTH - OFF T2', 'FontName','AvantGarde', 'FontSize',22,'FontWeight','bold', 'Color', [0 0 0]); 
set(gca,'FontName', 'AvantGarde', 'FontSize', 18, 'Box', 'off', 'TickDir','out', 'TickLength', [.02 .02],'XColor',[.1 .1 .1],'YColor',[.1 .1 .1], 'LineWidth', 1)
xlabel('seconds', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
ylabel('Firing Rate (spikes/second)', 'FontName','AvantGarde', 'FontSize', 18, 'Color', [0 0 0]);
set(gcf, 'Color', [1 1 1]);
legend('0331', '0409', '0605')
save_figure_pdf('/Users/sravi/Documents/My Documents/Grants/Sep 2013/Classf figures/', 'off-t2-same-pulse', gcf)





















































